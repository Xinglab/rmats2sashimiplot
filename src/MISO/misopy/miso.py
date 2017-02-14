# -*- mode: python; -*-
##
## Main MISO interface
##
import os
import csv
import time
import sys
import subprocess
from collections import defaultdict
import logging

import pysam

import misopy
import misopy.gff_utils as gff_utils
import misopy.as_events as as_events
import misopy.run_miso as run_miso
import misopy.misc_utils as misc_utils
import misopy.run_events_analysis as run_events
from misopy.parse_csv import *
from misopy.settings import Settings, load_settings
from misopy.settings import miso_path as miso_settings_path
import misopy.cluster_utils as cluster_utils

miso_path = os.path.dirname(os.path.abspath(__file__))
manual_url = "http://genes.mit.edu/burgelab/miso/docs/"


def get_main_logger(log_outdir,
                    level=logging.WARNING,
                    include_stdout=True):
    """
    Return logger object for main MISO thread.
    """
    logger_name = "miso_main"
    misc_utils.make_dir(log_outdir)
    logger = logging.getLogger(logger_name)
    formatter = \
        logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                          datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.root.setLevel(level)
    # Optionally add handler that streams all logs
    # to stdout
    if include_stdout:
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(level)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    # Write to main logger filename along
    # with time stamp
    logger_basename = "main.%s.log" %(misc_utils.get_timestamp())
    logger_fname = os.path.join(log_outdir, logger_basename)
    fh = logging.FileHandler(logger_fname)
    fh.setLevel(level)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger


def greeting(parser=None):
    print "MISO (Mixture of Isoforms model)"
    print "Probabilistic analysis of RNA-Seq data for detecting " \
          "differential isoforms"
    print "Use --help argument to view options.\n"
    if parser is not None:
        parser.print_help()
    

class GenesDispatcher:
    """
    Send MISO commands to cluster or locally
    using multi-cores.
    """
    def __init__(self, gff_dir, bam_filename,
                 output_dir, read_len, overhang_len,
                 main_logger,
                 settings_fname=None,
                 paired_end=None,
                 use_cluster=False,
                 chunk_jobs=200,
                 SGEarray=False,
                 sge_job_name="misojob",
                 gene_ids=None,
                 num_proc=None,
                 wait_on_jobs=True):
        self.main_logger = main_logger
        self.threads = {}
        self.gff_dir = gff_dir
        self.bam_filename = bam_filename
        # Check that the BAM filename exists and that it has an index
        if not os.path.isfile(self.bam_filename):
            self.main_logger.error("BAM file %s not found." %(self.bam_filename))
            sys.exit(1)
        self.bam_index_fname = "%s.bai" %(self.bam_filename)
        if not os.path.isfile(self.bam_index_fname):
            self.main_logger.warning("Expected BAM index file %s not found." \
                                %(self.bam_index_fname))
            self.main_logger.warning("Are you sure your BAM file is indexed?")
        self.output_dir = output_dir
        self.read_len = read_len
        # For now setting overhang to 1 always
        #self.overhang_len = overhang_len
        self.overhang_len = 1
        self.settings_fname = settings_fname
        self.paired_end = paired_end
        self.use_cluster = use_cluster
        self.chunk_jobs = chunk_jobs
        self.settings = Settings.get()
        self.cluster_cmd = Settings.get_cluster_command()
        self.sge_job_name = sge_job_name
        self.wait_on_jobs = wait_on_jobs
        # if chunk_jobs not given (i.e. set to False),
        # then set it to arbitrary value
        if not self.chunk_jobs:
            self.chunk_jobs = 200
        self.SGEarray = SGEarray
        self.num_processors = Settings.get_num_processors()
        if num_proc is not None:
            num_proc = int(num_proc)
            self.num_processors = num_proc
            self.main_logger.info("Using %d processors" %(num_proc))
        self.long_thresh = 50
        self.batch_logs_dir = \
            os.path.join(output_dir, "batch-logs")
        self.batch_genes_dir = \
            os.path.join(output_dir, "batch-genes")
        self.cluster_scripts_dir = \
            os.path.join(output_dir, "cluster_scripts")
        self.scripts_output_dir = \
            os.path.join(output_dir, "scripts_output")
        misc_utils.make_dir(self.batch_logs_dir)
        misc_utils.make_dir(self.batch_genes_dir)
        misc_utils.make_dir(self.cluster_scripts_dir)
        misc_utils.make_dir(self.scripts_output_dir)
        # First compile a set of genes that should be run on
        # and output them to file along with their indexed
        # filenames
        self.gene_ids_to_gff_index = \
            gff_utils.get_gene_ids_to_gff_index(gff_dir)
        # If we're given filtered gene IDs, use them
        if gene_ids is not None:
            self.gene_ids = gene_ids
        else:
            self.gene_ids = self.gene_ids_to_gff_index.keys()
        if len(self.gene_ids) == 0:
            self.main_logger.error("No genes to run on. Did you pass me the wrong path " \
                                   "to your index GFF directory? " \
                                   "Or perhaps your indexed GFF directory " \
                                   "is empty?")
            sys.exit(1)
        self.batch_filenames = self.output_batch_files()


    def output_batch_files(self):
        """
        Output a series of batch files containing
        gene IDs and their indexed GFF filenames.

        Return the batch filenames and their size.
        """
        batch_filenames = []
        chunk_jobs = self.chunk_jobs
        num_genes = len(self.gene_ids)
        num_chunks = max(1, int(round(num_genes / float(chunk_jobs))))
        if not self.use_cluster:
            # When not using cluster, use local multi-cores
            # using default number of processors
            num_chunks = self.num_processors
        gene_ids_batches = cluster_utils.chunk_list(self.gene_ids,
                                                    num_chunks)
        for batch_num, gene_ids_batch in enumerate(gene_ids_batches):
            batch_size = len(gene_ids_batch)
            batch_fname = os.path.join(self.batch_genes_dir,
                                       "batch-%d_genes.txt" %(batch_num))
            with open(batch_fname, "w") as batch_out:
                for gene_id in gene_ids_batch:
                    if gene_id not in self.gene_ids_to_gff_index:
                        # If gene is not found (perhaps because it had only a 'gene'
                        # entry in GFF, with no mRNA children), then skip it
                        print "Skipping: %s" %(gene_id)
                        continue
                    index_fname = self.gene_ids_to_gff_index[gene_id]
                    output_line = "%s\t%s\n" %(gene_id,
                                               index_fname)
                    batch_out.write(output_line)
            batch_filenames.append((batch_fname, batch_size))
        return batch_filenames


    def run(self, delay_constant=0.9):
        """
        Run batches either locally on multi-cores
        or using cluster.
        """
        batch_filenames = self.output_batch_files()
        # All MISO commands, each correspond to a batch,
        # and the number of jobs in each batch
        all_miso_cmds = []
        num_batches = len(batch_filenames)
        ##
        ## Prepare all the files necessary to run each batch
        ##
        print "Preparing to run %d batches of jobs..." %(num_batches)
        miso_run = os.path.join(miso_path, "run_miso.py")
        for batch_num, batch in enumerate(batch_filenames):
            batch_filename, batch_size = batch
            miso_cmd = \
              "python %s --compute-genes-from-file \"%s\" %s %s --read-len %d " \
                    %(miso_run,
                      batch_filename,
                      self.bam_filename,
                      self.output_dir,
                      self.read_len)
            # Add paired-end parameters and read len/overhang len
            if self.paired_end != None:
                # Run in paired-end mode
                frag_mean = float(self.paired_end[0])
                frag_sd = float(self.paired_end[1])
                miso_cmd += " --paired-end %.1f %.1f" %(frag_mean,
                                                        frag_sd)
            else:
                # Overhang len only used in single-end mode
                miso_cmd += " --overhang-len %d" %(self.overhang_len)
            # Add settings filename if given
            if self.settings_fname != None:
                miso_cmd += " --settings-filename %s" \
                    %(self.settings_fname)
            all_miso_cmds.append((miso_cmd, batch_size))
        ##
        ## Run all MISO commands for the batches
        ## either locally using multi-cores or on cluster
        ##
        # First handle special case of SGE cluster submission
        if self.use_cluster and self.SGEarray:
            print "Using SGEarray..."
            # Call SGE
            batch_argfile = os.path.join(self.cluster_scripts_dir,
                                         "run_args.txt")
            cluster_utils.run_SGEarray_cluster(all_miso_cmds,
                                               batch_argfile,
                                               self.output_dir,
                                               settings=self.settings_fname,
                                               job_name=self.sge_job_name,
                                               chunk=self.chunk_jobs)
            # End SGE case
            return
        # All cluster jobs 
        cluster_jobs = []
        for batch_num, cmd_info in enumerate(all_miso_cmds):
            miso_cmd, batch_size = cmd_info
            print "Running batch of %d genes.." %(batch_size)
            print "  - Executing: %s" %(miso_cmd)
            # Make a log file for the batch, where all the output
            # will be redirected
            time_str = time.strftime("%m-%d-%y_%H:%M:%S")
            batch_logfile = os.path.join(self.batch_logs_dir,
                                         "batch-%d-%s.log" %(batch_num,
                                                             time_str))
            cmd_to_run = "%s >> \"%s\";" %(miso_cmd, batch_logfile)
            if not self.use_cluster:
                # Run locally
                p = subprocess.Popen(cmd_to_run, shell=True)
                thread_id = "batch-%d" %(batch_num)
                print "  - Submitted thread %s" %(thread_id)
                self.threads[thread_id] = p
            else:
                # Run on cluster
                if batch_size >= self.long_thresh:
                    queue_type = "long"
                else:
                    queue_type = "short"
                # Run on cluster
                job_name = "gene_psi_batch_%d" %(batch_num)
                print "Submitting to cluster: %s" %(cmd_to_run)
                job_id = \
                    cluster_utils.run_on_cluster(cmd_to_run,
                                                 job_name,
                                                 self.output_dir,
                                                 queue_type=queue_type,
                                                 settings_fname=self.settings_fname)
                if job_id is not None:
                    cluster_jobs.append(job_id)
                time.sleep(delay_constant)
            # Extra delay constant
            time.sleep(delay_constant)                
        # If ran jobs on cluster, wait for them if there are any
        # to wait on.
        if self.wait_on_jobs:
            if self.use_cluster and (len(cluster_jobs) == 0):
                # If we're asked to use the cluster but the list
                # of cluster jobs is empty, it means we could not
                # find the IDs of the job from the submission
                # system. Report this to the user.
                self.main_logger.warning("Asked to wait on cluster jobs but cannot " \
                                         "parse their job IDs from the cluster submission " \
                                         "system.")
            # Try to wait on jobs no matter what; though if 'cluster_jobs'
            # is empty here, it will not wait
            cluster_utils.wait_on_jobs(cluster_jobs,
                                       self.cluster_cmd)
        else:
            if self.use_cluster:
                # If we're running in cluster mode and asked not
                # to wait for jobs, let user know
                self.main_logger.info("Not waiting on cluster jobs.")
        # If ran jobs locally, wait on them to finish
        # (this will do nothing if we submitted jobs to
        # cluster)
        self.wait_on_threads()
        

    def wait_on_threads(self):
        if self.use_cluster:
            # If ran jobs on cluster, nothing to wait for
            return
        threads_completed = {}
        num_threads = len(self.threads)
        if num_threads == 0:
            return
        print "Waiting on %d threads..." %(num_threads)
        t_start = time.time()
        for thread_name in self.threads:
            if thread_name in threads_completed:
                continue
            curr_thread = self.threads[thread_name]
            curr_thread.wait()
            if curr_thread.returncode != 0:
                self.main_logger.warning("Thread %s might have failed..." \
                                    %(thread_name))
            if curr_thread.returncode is None:
                self.main_logger.warning("Thread still going...")
            threads_completed[thread_name] = True
        t_end = time.time()
        duration = ((t_end - t_start) / 60.) / 60.
        self.main_logger.info("Threads completed in %.2f hours." \
                              %(duration))


def compute_all_genes_psi(gff_dir, bam_filename, read_len,
                          output_dir, main_logger,
                          use_cluster=False,
                          SGEarray=False,
                          chunk_jobs=800,
                          overhang_len=1,
                          paired_end=None,
                          settings_fname=None,
                          job_name="misojob",
                          num_proc=None,
                          prefilter=False,
                          wait_on_jobs=True):
    """
    Compute Psi values for genes using a GFF and a BAM filename.

    SGE functionality contributed by Michael Lovci.

    Options:
    - prefilter: if set to True, prefilter events by coverage.
      Uses bedtools to determine coverage of each event and remove
      events that do not meet the coverage criteria from the run.
    """
    print "Computing Psi values..." 
    print "  - GFF index: %s" %(gff_dir)
    print "  - BAM: %s" %(bam_filename)
    print "  - Read length: %d" %(read_len)
    print "  - Output directory: %s" %(output_dir)

    misc_utils.make_dir(output_dir)

    # Check GFF and BAM for various errors like headers mismatch
    run_events.check_gff_and_bam(gff_dir, bam_filename, main_logger,
                                 given_read_len=read_len)
    
    # Prefilter events that do not meet the coverage criteria
    # If filtering is on, only run on events that meet
    # the filter.
    all_gene_ids = None
    
    if prefilter:
        main_logger.info("Prefiltering on")
        if misc_utils.which("bedtools") is None:
            main_logger.error("Error: Cannot use bedtools. Bedtools is " \
                              "required for --prefilter option")
            sys.exit(1)
        filtered_gene_ids = run_events.get_ids_passing_filter(gff_dir,
                                                              bam_filename,
                                                              output_dir)
        # Prefiltering succeeded, so process only gene ids that
        # pass the filter
        if filtered_gene_ids != None:
            num_pass = len(filtered_gene_ids)
            all_gene_ids = filtered_gene_ids
            # If none of the events meet the read coverage filter
            # something must have gone wrong, e.g. mismatch
            # in chromosome headers between BAM and GFF
            if num_pass == 0:
                main_logger.error("None of the events in %s appear to meet the " \
                                  "read coverage filter. Check that your BAM headers " \
                                  "in %s match the GFF headers of indexed events." \
                                  %(gff_dir,
                                    bam_filename))
                sys.exit(1)
            main_logger.info("Total of %d events pass coverage filter." \
                             %(num_pass))

    ##
    ## Submit jobs either using cluster or locally
    ## using multi-cores.
    ##
    dispatcher = GenesDispatcher(gff_dir,
                                 bam_filename,
                                 output_dir,
                                 read_len,
                                 overhang_len,
                                 main_logger,
                                 settings_fname=settings_fname,
                                 paired_end=paired_end,
                                 use_cluster=use_cluster,
                                 chunk_jobs=chunk_jobs,
                                 sge_job_name=job_name,
                                 SGEarray=SGEarray,
                                 gene_ids=all_gene_ids,
                                 num_proc=num_proc,
                                 wait_on_jobs=wait_on_jobs)
    dispatcher.run()



def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--run", dest="compute_genes_psi",
                      nargs=2, default=None,
                      help="Compute Psi values for a given GFF annotation "
                      "of either whole mRNA isoforms or isoforms produced by "
                      "single alternative splicing events. Expects two "
                      "arguments: an indexed GFF directory with genes to "
                      "process, and a sorted, indexed BAM file (with "
                      "headers) to run on.")
    parser.add_option("--event-type", dest="event_type", nargs=1,
                      help="[OPTIONAL] Type of event (e.g. SE, RI, A3SS, ...)",
                      default=None)
    parser.add_option("--use-cluster", dest="use_cluster",
                      action="store_true", default=False,
                      help="Run events on cluster.")
    parser.add_option("--chunk-jobs", dest="chunk_jobs",
                      default=False, type="int",
                      help="Size (in number of events) of each job to chunk "
                      "events file into. Only applies when running on cluster.")
    parser.add_option("--no-filter-events", dest="no_filter_events",
                      action="store_true", default=False,
                      help="Do not filter events for computing Psi. "
                      "By default, MISO computes Psi only for events that "
                      "have a sufficient number of junction reads. "
                      "The default filter varies by event type.")
    parser.add_option("--settings-filename", dest="settings_filename",
                      default=os.path.join(miso_settings_path,
                                           "settings",
                                           "miso_settings.txt"),                    
                      help="Filename specifying MISO settings.")
    parser.add_option("--read-len", dest="read_len", default=None, type="int",
                      help="Length of sequenced reads.")
    parser.add_option("--paired-end", dest="paired_end", nargs=2, default=None,
                      help="Run in paired-end mode. Takes mean and "
                      "standard deviation of insert length distribution.")
    parser.add_option("--overhang-len", dest="overhang_len",
                      default=None, type="int",
                      help="Length of overhang constraints "
                      "imposed on junctions.")
    parser.add_option("--output-dir", dest="output_dir", default=None,
                      help="Directory for MISO output.")
    parser.add_option("--job-name", dest="job_name", nargs=1,
                      help="Name for jobs submitted to queue for SGE jobs. " \
                      "Default is misojob", default="misojob")
    parser.add_option("--SGEarray", dest="SGEarray",
                      action="store_true", default=False,
                      help="Use MISO on cluster with Sun Grid Engine. "
                      "To be used in conjunction with --use-cluster option.")
    parser.add_option("--prefilter", dest="prefilter", default=False,
                      action="store_true",
                      help="Prefilter events based on coverage. If given as " 
                      "argument, run will begin by mapping BAM reads to event "
                      "regions (using bedtools), and omit events that do not "
                      "meet coverage criteria from the run. By default, turned "
                      "off. Note that events that do not meet the coverage criteria "
                      "will not be processed regardless, but --prefilter simply "
                      "does this filtering step at the start of the run, potentially "
                      "saving computation time so that low coverage events will not "
                      "be processed or distributed to jobs if MISO is run on a "
                      "cluster. This options requires bedtools to be installed and "
                      "available on path.")
    parser.add_option("-p", dest="num_proc", default=None, nargs=1,
                      help="Number of processors to use. Only applies when running " \
                      "MISO on a single machine with multiple cores; does not apply " \
                      "to runs submitted to cluster with --use-cluster.")
    parser.add_option("--version", dest="version", default=False,
                      action="store_true",
                      help="Print MISO version.")
    parser.add_option("--no-wait", dest="no_wait", default=False,
                      action="store_true",
                      help="If passed in, do not wait on cluster jobs after " \
                      "they are submitted. By default, wait.")
    ##
    ## Gene utilities
    ##
    parser.add_option("--view-gene", dest="view_gene",
                      nargs=1, default=None,
                      help="View the contents of a gene/event that has "
                      "been indexed. Takes as input an "
                      "indexed (.pickle) filename.")
    (options, args) = parser.parse_args()

    greeting()

    if options.version:
        print "MISO version %s\n" %(misopy.__version__)

    ##
    ## Load the settings file 
    ##
    if not os.path.isdir(miso_settings_path):
        print "Error: %s is not a directory containing a default MISO " \
              "settings filename. Please specify a settings filename " \
              "using --settings-filename."
        return
    
    settings_filename = \
        os.path.abspath(os.path.expanduser(options.settings_filename))
    Settings.load(settings_filename)
    
    if (not options.use_cluster) and options.chunk_jobs:
        print "Error: Chunking jobs only applies when using " \
              "the --use-cluster option to run MISO on cluster."
        sys.exit(1)
    if (not options.use_cluster) and options.SGEarray:
        print "Error: SGEarray implies that you are using an SGE cluster," \
              "please run again with --use-cluster option enabled."
        sys.exit(1)

    ##
    ## Quantitation using BAM for all genes
    ##
    if options.compute_genes_psi != None:
        # GFF filename with genes to process
        gff_filename = \
            os.path.abspath(os.path.expanduser(options.compute_genes_psi[0]))

        # BAM filename with reads
        bam_filename = \
            os.path.abspath(os.path.expanduser(options.compute_genes_psi[1]))

        if options.output_dir == None:
            print "Error: need --output-dir to compute Psi values."
            sys.exit(1)

        # Output directory to use
        output_dir = os.path.abspath(os.path.expanduser(options.output_dir))

        ##
        ## Load the main logging object
        ##
        logs_output_dir = os.path.join(output_dir, "logs")
        main_logger = get_main_logger(logs_output_dir)

        if options.read_len == None:
            main_logger.error("need --read-len to compute Psi values.")
            sys.exit(1)

        overhang_len = 1

        if options.paired_end != None and options.overhang_len != None:
            main_logger.warning("cannot use --overhang-len in paired-end mode.\n" \
                                "Using overhang = 1")
        if options.overhang_len != None:
            overhang_len = options.overhang_len

        # Whether to wait on cluster jobs or not
        wait_on_jobs = not options.no_wait
        compute_all_genes_psi(gff_filename, bam_filename,
                              options.read_len, output_dir,
                              main_logger,
                              overhang_len=overhang_len,
                              use_cluster=options.use_cluster,
                              SGEarray=options.SGEarray,
                              job_name=options.job_name,
                              chunk_jobs=options.chunk_jobs,
                              paired_end=options.paired_end,
                              settings_fname=settings_filename,
                              prefilter=options.prefilter,
                              num_proc=options.num_proc,
                              wait_on_jobs=wait_on_jobs)

    if options.view_gene != None:
        indexed_gene_filename = \
            os.path.abspath(os.path.expanduser(options.view_gene))
        print "Viewing genes in %s" %(indexed_gene_filename)
        gff_genes = gff_utils.load_indexed_gff_file(indexed_gene_filename)

        if gff_genes == None:
            print "No genes."
            sys.exit(1)

        for gene_id, gene_info in gff_genes.iteritems():
            print "Gene %s" %(gene_id)
            gene_obj = gene_info['gene_object']
            print " - Gene object: ", gene_obj
            print "=="
            print "Isoforms: "
            for isoform in gene_obj.isoforms:
                print " - ", isoform
            print "=="
            print "mRNA IDs: "
            for mRNA_id in gene_info['hierarchy'][gene_id]['mRNAs']:
                print "%s" %(mRNA_id)
            print "=="    
            print "Exons: "
            for exon in gene_obj.parts:
                print " - ", exon
    
if __name__ == "__main__":
    main()
