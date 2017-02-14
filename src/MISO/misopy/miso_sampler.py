##
## MCMC sampler for Mixture-of-Isoforms (MISO) model
##
## Yarden Katz <yarden@mit.edu>
##
## The sampler uses a Metropolis-Hastings sampling scheme, combined with
## a Gibbs sampling step.
##
import scipy
import misopy
from misopy.reads_utils import count_aligned_reads, \
                               count_isoform_assignments
from misopy.read_simulator import simulate_reads, print_reads_summary, \
                                  read_counts_to_read_list, \
                                  get_reads_summary
import misopy.hypothesis_test as ht
from misopy.Gene import Gene, Exon
from misopy.py2c_gene import *

# C MISO interface
import pysplicing

from scipy import *
from numpy import *
import cPickle as pickle
from scipy.stats import mode
import math
import time
from numpy import numarray
import os
import sys
from collections import defaultdict
import glob
import logging
import logging.handlers

loggers = {}

def get_logger(logger_name, log_outdir,
               level=logging.WARNING,
               include_stdout=True):
    """
    Return a logging object.
    """
    global loggers
    # Avoid race-conditions
    try:
        os.makedirs(log_outdir)
    except OSError:
        pass
    if loggers.get(logger_name):
        return loggers.get(logger_name)
    logger = logging.getLogger(logger_name)
    formatter = \
        logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                          datefmt='%m/%d/%Y %I:%M:%S %p')
    # Do not log to file
    #if log_outdir is not None:
    #    log_filename = os.path.join(log_outdir, "%s.log" %(logger_name))
    #    fh = logging.FileHandler(log_filename)
    #    fh.setLevel(level)
    #    fh.setFormatter(formatter)
    #    logger.addHandler(fh)
    logging.root.setLevel(level)
    # Optionally add handler that streams all logs
    # to stdout
    if include_stdout:
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(level)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    logger.info("Created logger %s" %(logger_name))
    loggers.update({logger_name: logger})
    return logger


##
## Helper statistics/linear algebra functions
##
def set_diag(a, v):
    for i, elt in enumerate(a):
        a[i, i] = v
    return a

def maxi(l):
    m = max(l)
    for i, v in enumerate(l):
        if m == v:
            return i

def mini(l):
    m = min(l)
    for i, v in enumerate(l):
        if m == v:
            return i


def exp_logsumexp(a):
    return exp(a - logsumexp(a))


def vect_logsumexp(a, axis=None):
    if axis is None:
        # Use the scipy.maxentropy version.
        return logsumexp(a)
    a = asarray(a)
    shp = list(a.shape)
    shp[axis] = 1
    a_max = a.max(axis=axis)
    s = log(exp(a - a_max.reshape(shp)).sum(axis=axis))
    lse  = a_max + s
    return lse


def print_assignment_summary(assignments):
    counts = defaultdict(int)
    for a in assignments:
        counts[a] += 1
    for k, v in counts.iteritems():
        print "Total of %d in isoform %d" %(v, k)


def float_array_to_str(array_of_floats):
    """
    Convert a float numpy array to a string for printing purposes.
    """
    str_float_array = '[' + ' '.join(['%.3f' %(val) for val in array_of_floats]) + ']'
    return str_float_array


def get_paired_end_sampler_params(num_isoforms,
                                  mean_frag_len,
                                  frag_variance,
                                  read_len,
                                  overhang_len=1):
    """
    Return parameters for MISO sampler, in paired-end mode.
    """
    hyperparameters = ones(num_isoforms)
    proposal_diag = 0.05
    sigma = set_diag(zeros([num_isoforms-1, num_isoforms-1]),
                     proposal_diag)
    sampler_params = {'read_len': read_len,
                      'overhang_len': overhang_len,
                      'uniform_proposal': False,
                      'sigma_proposal': sigma,
                      'mean_frag_len': mean_frag_len,
                      'frag_variance': frag_variance}
    return sampler_params


def get_single_end_sampler_params(num_isoforms,
                                  read_len,
                                  overhang_len=1):
    """
    Return parameters for MISO sampler, in single-end mode.
    """
    hyperparameters = ones(num_isoforms)
    proposal_diag = 0.05
    sigma = set_diag(zeros([num_isoforms-1, num_isoforms-1]),
                     proposal_diag)
    sampler_params = {'read_len': read_len,
                      'overhang_len': overhang_len,
                      'uniform_proposal': False,
                      'sigma_proposal': sigma}
    return sampler_params


class MISOSampler:
    def __init__(self, params,
                 paired_end=False,
                 log_dir=None):
        """
        Make a sampler with the given parameters.
        """
        self.params = params
        self.paired_end = paired_end
        # set default fragment length distribution parameters
        if self.paired_end:
            if ((not 'mean_frag_len' in self.params) or \
                (not 'frag_variance' in self.params)):
                raise Exception, "Must set mean_frag_len and frag_variance when " \
                      "running in sampler on paired-end data."
            self.mean_frag_len = self.params['mean_frag_len']
            self.frag_variance = self.params['frag_variance']

        if log_dir != None:
            self.log_dir = os.path.abspath(os.path.expanduser(log_dir))
            self.log_dir = os.path.join(log_dir, 'logs')
            # Avoid race-conditions
            try:
                os.makedirs(self.log_dir)
            except OSError:
                pass
        self.miso_logger = get_logger('miso_logger', self.log_dir)
        self.miso_logger.info("Instantiated sampler.")


    def run_sampler(self, num_iters, reads, gene, hyperparameters, params,
                    output_file,
                    num_chains=6,
                    burn_in=1000,
                    lag=2,
                    prior_params=None,
                    # By default, use sampler with read classes (collapsed)
                    # to get speed boost for single-end reads
                    # (To revert to old reassigning sampler, use
                    # pysplicing.MISO_ALGO_REASSIGN)
                    algorithm=pysplicing.MISO_ALGO_CLASSES,
                    start_cond=pysplicing.MISO_START_AUTO,
                    stop_cond=pysplicing.MISO_STOP_FIXEDNO,
                    verbose=True):
        """
        Fast version of MISO MCMC sampler.

        Calls C version and returns results.
        """
        num_isoforms = len(gene.isoforms)
        self.num_isoforms = num_isoforms

        if prior_params == None:
            prior_params = (1.0,) * num_isoforms

        read_positions = reads[0]
        read_cigars = reads[1]

        self.num_reads = len(read_positions)

        if self.num_reads == 0:
            print "No reads for gene: %s" %(gene.label)
            return

        output_file = output_file + ".miso"
        # If output filename exists, don't run sampler
        if os.path.isfile(os.path.normpath(output_file)):
            print "Output filename %s exists, not running MISO." \
                %(output_file)
            return None

        self.params['iters'] = num_iters
        self.params['burn_in'] = burn_in
        self.params['lag'] = lag

        # Define local variables related to reads and overhang
        self.overhang_len = self.params['overhang_len']
        self.read_len = self.params['read_len']

        t1 = 0
        t2 = 0
        if verbose:
            t1 = time.time()
        #self.miso_logger.info("Running sampler...")
        #self.miso_logger.info("  - num_iters: " + str(num_iters))
        #self.miso_logger.info("  - burn-in: " + str(burn_in))
        #self.miso_logger.info("  - lag: " + str(lag))
        #self.miso_logger.info("  - paired-end? " + str(self.paired_end))
        #self.miso_logger.info("  - gene: " + str(gene))
        rejected_proposals = 0
        accepted_proposals = 0
        psi_vectors = []
        all_psi_proposals = []

        if params['uniform_proposal']:
            self.miso_logger.debug("UNIFORM independent proposal being used.")
            proposal_type = "unif"
        else:
            self.miso_logger.debug("Non-uniform proposal being used.")
            self.miso_logger.debug("  - sigma_proposal: " + str(params['sigma_proposal']))
            proposal_type = "drift"
        init_psi = ones(num_isoforms)/float(num_isoforms)
        # Do not process genes with one isoform
        if num_isoforms == 1:
            one_iso_msg = "Gene %s has only one isoform; skipping..." \
                          %(gene.label)
            self.miso_logger.warning(one_iso_msg)
            return

        # Convert Python Gene object to C
        c_gene = py2c_gene(gene)

        ##
        ## Run C MISO
        ##
        read_positions = tuple([r+1 for r in read_positions])
        if self.paired_end:
            # Number of standard deviations in insert length
            # distribution to consider when assigning reads
            # to isoforms
            num_sds = 4L

            # Run paired-end
            miso_results = pysplicing.MISOPaired(c_gene, 0L,
                                                 read_positions,
                                                 read_cigars,
                                                 long(self.read_len),
                                                 float(self.mean_frag_len),
                                                 float(self.frag_variance),
                                                 float(num_sds),
                                                 long(num_iters),
                                                 long(burn_in),
                                                 long(lag),
                                                 prior_params,
                                                 long(self.overhang_len),
                                                 long(num_chains),
                                                 start_cond,
                                                 stop_cond)
        else:
            # Run single-end
            miso_results = pysplicing.MISO(c_gene,
                                           0L,
                                           read_positions,
                                           read_cigars,
                                           long(self.read_len),
                                           long(num_iters),
                                           long(burn_in),
                                           long(lag),
                                           prior_params,
                                           long(self.overhang_len),
                                           long(num_chains),
                                           start_cond,
                                           stop_cond,
                                           pysplicing.MISO_ALGO_REASSIGN)
#                                           algorithm)

        # Psi samples
        psi_vectors = transpose(array(miso_results[0]))

        # Log scores of accepted samples
        kept_log_scores = transpose(array(miso_results[1]))

        # Read classes
        read_classes = miso_results[2]

        # Read class statistics
        read_class_data = miso_results[3]

        # Assignments of reads to isoforms
        assignments = miso_results[4]

        # Statistics and parameters about sampler run
        run_stats = miso_results[5]

        # Assignments of reads to classes.
        # read_classes[n] represents the read class that has
        # read_assignments[n]-many reads.
        reads_data = (read_classes, read_class_data)

        assignments = array(assignments)

        # Skip events where all reads are incompatible with the annotation;
        # do not output a file for those.
        if all(assignments == -1):
            print "All reads incompatible with annotation, skipping..."
            return

        accepted_proposals = run_stats[4]
        rejected_proposals = run_stats[5]

        percent_acceptance = (float(accepted_proposals)/(accepted_proposals + \
                                                         rejected_proposals)) * 100
        #self.miso_logger.info("Percent acceptance (including burn-in): %.4f" %(percent_acceptance))
        #self.miso_logger.info("Number of iterations recorded: %d" %(len(psi_vectors)))

        # Write MISO output to file
        print "Outputting samples to: %s..." %(output_file)
        self.miso_logger.info("Outputting samples to: %s" %(output_file))
        self.output_miso_results(output_file, gene, reads_data, assignments,
                                 psi_vectors, kept_log_scores, num_iters,
                                 burn_in, lag, percent_acceptance,
                                 proposal_type)
        if verbose:
            t2 = time.time()
            print "Event took %.2f seconds" %(t2 - t1)


    def output_miso_results(self, output_file, gene, reads_data, assignments,
                            psi_vectors, kept_log_scores, num_iters, burn_in,
                            lag, percent_acceptance, proposal_type):
        """
        Output results of MISO to a file.
        """
        output = open(output_file, 'w')

        # Get a string representation of the isoforms - use '_'
        # in the delimiter regardless
        iso_delim = '_'
        if type(gene.isoforms[0].desc) == list:
            str_isoforms = '[' + ",".join(["\'" + iso_delim.join(iso.desc) + "\'" \
                                           for iso in gene.isoforms]) + ']'
        else:
            str_isoforms = '[' + ",".join(["\'" + iso.desc + "\'" \
                                           for iso in gene.isoforms]) + ']'

        num_isoforms = len(gene.isoforms)

        # And of the exon lengths
        exon_lens = ",".join(["(\'%s\',%d)" %(p.label, p.len) \
                              for p in gene.parts])

        ## Compile header with information about isoforms and internal parameters used
        ## by the sampler, and also information about read counts and number of
        ## reads assigned to each isoform.

        read_classes, read_class_counts = reads_data
        read_counts_list = []

        for class_num, class_type in enumerate(read_classes):
            class_counts = read_class_counts[class_num]

            # Get the read class type in string format
            class_str = str(tuple([int(c) for c in class_type])).replace(" ", "")

            # Get the read class counts in string format
            class_counts_str = "%s" %(int(read_class_counts[class_num]))

            # Put class and counts together
            curr_str = "%s:%s" %(class_str,
                                 class_counts_str)
            read_counts_list.append(curr_str)

        # Get a summary of the raw read counts supporting each isoform
        read_counts_str = ",".join(read_counts_list)

        assigned_counts = count_isoform_assignments(assignments)

        # Get number of reads assigned to each isoform
        assigned_counts_str = ",".join(["%d:%d" %(c[0], c[1]) \
                                        for c in assigned_counts])

        # coordinates where mRNAs start
        mRNA_starts = []
        mRNA_ends = []
        for iso in gene.isoforms:
            mRNA_starts.append(iso.genomic_start)
            mRNA_ends.append(iso.genomic_end)
        mRNA_start_coords = ",".join([str(start) for start in mRNA_starts])
        mRNA_end_coords = ",".join([str(end) for end in mRNA_ends])
        chrom = gene.chrom
        if chrom == None:
            chrom = "NA"
        strand = gene.strand
        if strand == None:
            strand = "NA"
        header = "#isoforms=%s\texon_lens=%s\titers=%d\tburn_in=%d\tlag=%d\t" \
                 "percent_accept=%.2f\tproposal_type=%s\t" \
                 "counts=%s\tassigned_counts=%s\tchrom=%s\tstrand=%s\tmRNA_starts=%s\tmRNA_ends=%s\n" \
                 %(str_isoforms, exon_lens, num_iters, burn_in, lag,
                   percent_acceptance, proposal_type, read_counts_str,
                   assigned_counts_str,
                   # Fields related to gene/event
                   chrom,
                   strand,
                   mRNA_start_coords,
                   mRNA_end_coords)
        output.write(header)

        # Output samples and their associated log scores, as well as read counts
        results_fields = ["sampled_psi", "log_score"]
        results_header = "%s\n" %("\t".join(results_fields))
        output.write(results_header)
        for psi_sample, curr_log_score in zip(psi_vectors, kept_log_scores):
            psi_sample_str = ",".join(["%.4f" %(psi) for psi in psi_sample])
            output_line = "%s\t%.2f\n" %(psi_sample_str, curr_log_score)
            output.write(output_line)
        output.close()
        print "Completed outputting."
#        return [percent_acceptance, array(psi_vectors), array(kept_log_scores)]

def run_sampler_on_event(gene, ni, ne, nb, read_len, overhang_len, num_iters,
                         output_dir, confidence_level=.95):
    """
    Run sampler on a two-isoform gene event.
    """
    print "Running sampler on a two-isoform event..."
    print "  - Gene label: ", gene.label, gene
    print "  - NI, NE, NB: %d, %d, %d" %(ni, ne, nb)
    print "Using default sampler parameters."
    if gene.chrom != None:
        # Index output by chromosome
        print "Indexing by chromosome..."
        output_dir = os.path.join(output_dir, gene.chrom)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    output_filename = os.path.join(output_dir, gene.label)

    samples = []
    cred_interval = []

    num_isoforms = len(gene.isoforms)
    burn_in = 500
    lag = 10
    hyperparameters = ones(num_isoforms)
    proposal_diag = 0.05
    sigma = set_diag(zeros([num_isoforms-1, num_isoforms-1]),
                     proposal_diag)
    sampler_params = {'read_len': read_len,
                      'overhang_len': overhang_len,
                      'uniform_proposal': False,
                      'sigma_proposal': sigma}
    sampler = MISOSampler(sampler_params, log_dir=output_dir)
    reads = read_counts_to_read_list(ni, ne, nb)
    t1 = time.time()
    sampler_results = sampler.run_sampler(num_iters, reads, gene, hyperparameters,
                                          sampler_params, output_filename, burn_in=burn_in,
                                          lag=lag)
    if not sampler_results:
        return (samples, cred_interval)
    samples = sampler_results[1]
    # Compute credible intervals
    cred_interval = ht.compute_credible_intervals(samples, confidence_level=confidence_level)
    t2 = time.time()
    print "  - Sampler run took %s seconds." %(str(t2-t1))
    # return samples and credible intervals
    return (samples, cred_interval)


def profile_miso():
    from Gene import make_gene
    gene = make_gene([150, 100, 150], [[1, 2, 3], [1, 3]])
    read_len = 36
    overhang_len = 4
    output_dir = "profiler-test"
    for x in range(10):
        print "x = %d" %(x)
        a, b = run_sampler_on_event(gene, 500, 50, 40, read_len, overhang_len,
                                    10000, output_dir)



def main():
    return
    # import cProfile as profile
    # import pstats
    # output_file = "profile"
    # profile.run('profile_miso()', output_file)
    # p = pstats.Stats(output_file)
    # print "name: "
    # print p.sort_stats('name')
    # print "all stats: "
    # p.print_stats()
    # print "cumulative (top 10): "
    # p.sort_stats('cumulative').print_stats(20)

if __name__ == '__main__':
    main()
