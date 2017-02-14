#!/usr/bin/env python
##
## Interface for running MISO locally or on cluster
##

import os
import sys
import time
import glob

import misopy
from misopy.settings import Settings
from misopy.settings import miso_path as miso_settings_path
import misopy.hypothesis_test as ht
import misopy.as_events as as_events
import misopy.cluster_utils as cluster_utils
import misopy.sam_utils as sam_utils
import misopy.miso_sampler as miso
import misopy.Gene as gene_utils
import misopy.gff_utils as gff_utils
import misopy.misc_utils as misc_utils

from misopy.parse_csv import *
from misopy.samples_utils import *

import numpy as np
np.seterr(all='ignore')

miso_path = os.path.dirname(os.path.abspath(__file__))

##
## General interface
##
def compute_gene_psi(gene_ids, gff_index_filename, bam_filename,
                     output_dir, read_len, overhang_len,
                     paired_end=None,
                     event_type=None,
                     verbose=True):
    """
    Run Psi at the Gene-level (for multi-isoform inference.)

    Arguments:

    - Set of gene IDs corresponding to gene IDs from the GFF
    - Indexed GFF filename describing the genes
    - BAM filename with the reads (must be sorted and indexed)
    - Output directory
    - Optional: Run in paired-end mode. Gives mean and standard deviation
      of fragment length distribution.
    """
    misc_utils.make_dir(output_dir)
        
    if not os.path.exists(gff_index_filename):
        print "Error: No GFF %s" %(gff_index_filename)
        return
    
    num_genes = len(gene_ids)
    
    print "Computing Psi for %d genes..." %(num_genes)
    print "  - " + ", ".join(gene_ids)
    print "  - GFF filename: %s" %(gff_index_filename)
    print "  - BAM: %s" %(bam_filename)
    print "  - Outputting to: %s" %(output_dir)

    if paired_end:
        print "  - Paired-end mode: ", paired_end

    settings = Settings.get()
    settings_params = Settings.get_sampler_params()
    burn_in = settings_params["burn_in"]
    lag = settings_params["lag"]
    num_iters = settings_params["num_iters"]
    num_chains = settings_params["num_chains"]

    min_event_reads = Settings.get_min_event_reads()
    strand_rule = Settings.get_strand_param()

    mean_frag_len = None
    frag_variance = None

    if paired_end:
        mean_frag_len = int(paired_end[0])
        frag_variance = power(int(paired_end[1]), 2)

    # Load the genes from the GFF
    gff_genes = gff_utils.load_indexed_gff_file(gff_index_filename)
    
    # If given a template for the SAM file, use it
    template = None

    if settings and "sam_template" in settings:
        template = settings["sam_template"]

    if "filter_reads" not in settings:
        filter_reads = True
    else:
        filter_reads = settings["filter_reads"]
        
    # Load the BAM file upfront
    bamfile = sam_utils.load_bam_reads(bam_filename,
                                       template=template)
    # Check if we're in compressed mode
    compressed_mode = misc_utils.is_compressed_index(gff_index_filename)
    
    for gene_id, gene_info in gff_genes.iteritems():
        lookup_id = gene_id
        # Skip genes that we were not asked to run on
        if lookup_id not in gene_ids:
            continue
        gene_obj = gene_info['gene_object']
        gene_hierarchy = gene_info['hierarchy']

        # Sanity check: if the isoforms are all shorter than the read,
        # skip the event
        if all(map(lambda l: l < read_len, gene_obj.iso_lens)):
            print "All isoforms of %s shorter than %d, so skipping" \
                  %(gene_id, read_len)
            continue
        
        # Find the most inclusive transcription start and end sites
        # for each gene
        tx_start, tx_end = \
            gff_utils.get_inclusive_txn_bounds(gene_info['hierarchy'][gene_id])

        # Fetch reads aligning to the gene boundaries
        gene_reads = \
            sam_utils.fetch_bam_reads_in_gene(bamfile,
                                              gene_obj.chrom,
                                              tx_start,
                                              tx_end,
                                              gene_obj)
        # Parse reads: checking strandedness and pairing
        # reads in case of paired-end data
        reads, num_raw_reads = \
            sam_utils.sam_parse_reads(gene_reads,
                                      paired_end=paired_end,
                                      strand_rule=strand_rule,
                                      target_strand=gene_obj.strand,
                                      given_read_len=read_len)
        # Skip gene if none of the reads align to gene boundaries
        if filter_reads:
            if num_raw_reads < min_event_reads:
                print "Only %d reads in gene, skipping (needed >= %d reads)" \
                      %(num_raw_reads,
                        min_event_reads)
                continue
            else:
                print "%d raw reads in event" %(num_raw_reads)

        num_isoforms = len(gene_obj.isoforms)
        hyperparameters = ones(num_isoforms)

        ##
        ## Run the sampler
        ##
        # Create the sampler with the right parameters depending on whether
        # this is a paired-end or single-end data set.
        if paired_end:
            # Sampler parameters for paired-end mode
            sampler_params = \
                miso.get_paired_end_sampler_params(num_isoforms,
                                                   mean_frag_len,
                                                   frag_variance,
                                                   read_len,
                                                   overhang_len=overhang_len)
            sampler = miso.MISOSampler(sampler_params,
                                       paired_end=True,
                                       log_dir=output_dir)

        else:
            # Sampler parameters for single-end mode
            sampler_params = miso.get_single_end_sampler_params(num_isoforms,
                                                                read_len,
                                                                overhang_len)
            sampler = miso.MISOSampler(sampler_params,
                                       paired_end=False,
                                       log_dir=output_dir)

        # Make directory for chromosome -- if given an event type, put
        # the gene in the event type directory
        if event_type != None:
            chrom_dir = os.path.join(output_dir, event_type, gene_obj.chrom)
        else:
            chrom_dir = os.path.join(output_dir, gene_obj.chrom)

        try:
            os.makedirs(chrom_dir)
        except OSError:
            pass

        # Pick .miso output filename based on the pickle filename
        miso_basename = os.path.basename(gff_index_filename)
        if not miso_basename.endswith(".pickle"):
            print "Error: Invalid index file %s" %(gff_index_filename)
            sys.exit(1)
        miso_basename = miso_basename.replace(".pickle", "")
        output_filename = os.path.join(chrom_dir, "%s" %(miso_basename))
        sampler.run_sampler(num_iters, reads, gene_obj, hyperparameters,
                            sampler_params, output_filename,
                            num_chains=num_chains,
                            burn_in=burn_in,
                            lag=lag)


def run_compute_genes_from_file(options):
    """
    Run on a set of genes/events described a file.

    File is two-column, tab-delimited where first column
    is the name of the event/gene (ID= from GFF) and the
    second column is the path to the indexed GFF event
    corresponding to the event/gene.
    """
    if options.read_len == None:
        print "Error: must provide --read-len."
        sys.exit(1)

    overhang_len = 1
    if options.overhang_len != None:
        overhang_len = options.overhang_len

    paired_end = None
    # Parse arguments from options
    genes_filename = \
        os.path.abspath(os.path.expanduser(options.compute_genes_from_file[0]))
    bam_filename = \
        os.path.abspath(os.path.expanduser(options.compute_genes_from_file[1]))
    output_dir = \
        os.path.abspath(os.path.expanduser(options.compute_genes_from_file[2]))
    print "Computing Psi for genes from file..."
    print "  - Input file: %s" %(genes_filename)
    if options.paired_end != None:
        paired_end = float(options.paired_end[0]), \
                     float(options.paired_end[1])
        print "  - Paired-end mode"
    # Check that the events filename exists
    if not os.path.isfile(genes_filename):
        print "Error: %s filename does not exist." %(genes_filename)
        sys.exit(1)
    if not os.path.isfile(bam_filename):
        print "Error: BAM filename %s does not exist." %(bam_filename)
        sys.exit(1)
    # Load the events and their indexed GFF paths
    num_genes = 0
    with open(genes_filename) as genes_in:
        for line in genes_in:
            gene_id, gff_filename = line.strip().split("\t")
            if not os.path.isfile(gff_filename):
                print "Error: %s does not exist." %(gff_filename)
                sys.exit(1)
            compute_gene_psi([gene_id], gff_filename, bam_filename,
                             output_dir, options.read_len, overhang_len,
                             paired_end=paired_end,
                             event_type=options.event_type)
            num_genes += 1
    print "Processed %d genes" %(num_genes)
            

def run_compute_gene_psi(options):
    """
    Parse options and run compute_genes_psi.
    """
    if options.read_len == None:
        print "Error: must provide --read-len."
        sys.exit(1)

    overhang_len = 1
    if options.overhang_len != None:
        overhang_len = options.overhang_len

    paired_end = None
    if options.paired_end != None:
        paired_end = float(options.paired_end[0]), \
                     float(options.paired_end[1])

    # Genes to run on from GFF
    gene_ids = options.compute_gene_psi[0].split(",")

    # GFF filename describing genes
    gff_filename = \
        os.path.abspath(os.path.expanduser(options.compute_gene_psi[1]))

    # BAM filename with reads
    bam_filename = \
        os.path.abspath(os.path.expanduser(options.compute_gene_psi[2]))

    # Output directory
    output_dir = \
        os.path.abspath(os.path.expanduser(options.compute_gene_psi[3]))

    compute_gene_psi(gene_ids, gff_filename, bam_filename, output_dir,
                     options.read_len, overhang_len,
                     paired_end=paired_end,
                     event_type=options.event_type)

        
def greeting(parser=None):
    print "MISO (Mixture of Isoforms model)"
    print "Probabilistic analysis of RNA-Seq data to detect " \
          "differential isoforms"
    print "Use --help argument to view options.\n"
    if parser is not None:
        parser.print_help()
    
    
def main():
    from optparse import OptionParser
    parser = OptionParser()

    ##
    ## Main options
    ##
    parser.add_option("--compute-gene-psi", dest="compute_gene_psi",
                      nargs=4, default=None,
                      help="Compute Psi using for a given multi-isoform gene. "
                      "Expects four arguments: the first is a gene ID or set "
                      "of comma-separated (no spaces) gene IDs, "
                      "the second is a GFF indexed file with the gene "
                      "information, the third is a sorted and "
                      "indexed BAM file with reads aligned to the gene, "
                      "and the fourth is an output directory.")
    parser.add_option("--paired-end", dest="paired_end",
                      nargs=2, default=None,
                      help="Run in paired-end mode.  Takes a mean and standard "
                      "deviation for the fragment length distribution (assumed "
                      "to have discretized normal form.)")
    parser.add_option("--compute-genes-from-file", dest="compute_genes_from_file",
                      nargs=3, default=None,
                      help="Runs on a set of genes from a file. Takes as input: "
                      "(1) a two-column tab-delimited file, where column 1 is the "
                      "event ID (ID field from GFF) and the second column is "
                      "the path to the indexed GFF file for that event. "
                      "MISO will run on all the events described in the file, "
                      "(2) a sorted, indexed BAM file to run on, and (3) a "
                      "directory to output results to.")
    
    ##
    ## Psi utilities
    ##
    parser.add_option("--compare-samples", dest="samples_to_compare",
                      nargs=3, default=None,
		      help="Compute comparison statistics between the two "
                      "given samples. Expects three directories: the first is "
                      "sample1's MISO output, the second is sample2's MISO "
                      "output, and the third is the directory where "
		      "results of the sample comparison will be outputted.")
    parser.add_option("--comparison-labels", dest="comparison_labels",
                      nargs=2, default=None,
                      help="Use these labels for the sample comparison "
                      "made by --compare-samples. "
                      "Takes two arguments: the label for sample 1 "
                      "and the label for sample 2, where sample 1 and "
                      "sample 2 correspond to the order of samples given "
                      "to --compare-samples.")
    parser.add_option("--summarize-samples", dest="summarize_samples",
                      nargs=2, default=None,
		      help="Compute summary statistics of the given set "
                      "of samples. Expects a directory with MISO output "
                      "and a directory to output summary file to.")
    parser.add_option("--summary-label", dest="summary_label",
                      nargs=1, default=None,
                      help="Label for MISO summary file. If not given, "
                      "uses basename of MISO output directory.")
    parser.add_option("--use-cluster", action="store_true",
                      dest="use_cluster", default=False)
    parser.add_option("--chunk-jobs", dest="chunk_jobs",
                      default=False, type="int",
		      help="Size (in number of events) of each job to "
                      "chunk events file into. Only applies when "
                      "running on cluster.")
    parser.add_option("--settings-filename", dest="settings_filename",
                      default=os.path.join(miso_settings_path,
                                           "settings",
                                           "miso_settings.txt"),
                      help="Filename specifying MISO settings.")
    parser.add_option("--read-len", dest="read_len", type="int",
                      default=None)
    parser.add_option("--overhang-len", dest="overhang_len", type="int",
                      default=None)
    parser.add_option("--event-type", dest="event_type", default=None,
		      help="Event type of two-isoform "
                      "events (e.g. 'SE', 'RI', 'A3SS', ...)")    
    parser.add_option("--use-compressed", dest="use_compressed",
                      nargs=1, default=None,
                      help="Use compressed event IDs. Takes as input a "
                      "genes_to_filenames.shelve file produced by the "
                      "index_gff script.")
    ##
    ## Gene utilities
    ##
    parser.add_option("--view-gene", dest="view_gene",
                      nargs=1, default=None,
                      help="View the contents of a gene/event that has "
                      "been indexed. Takes as input an "
                      "indexed (.pickle) filename.")
    (options, args) = parser.parse_args()

    if options.compute_gene_psi is None:
        greeting()

    ##
    ## Load the settings file 
    ##
    Settings.load(os.path.expanduser(options.settings_filename))

    use_compressed = None
    if options.use_compressed is not None:
        use_compressed = \
            os.path.abspath(os.path.expanduser(options.use_compressed))
        if not os.path.exists(use_compressed):
            print "Error: mapping filename from event IDs to compressed IDs %s " \
                  "is not found." %(use_compressed)
            sys.exit(1)
        else:
            print "Compression being used."
            
    if options.samples_to_compare is not None:
        sample1_dirname = os.path.abspath(options.samples_to_compare[0])
	sample2_dirname = os.path.abspath(options.samples_to_compare[1])
	output_dirname = os.path.abspath(options.samples_to_compare[2])
	if not os.path.isdir(output_dirname):
            print "Making comparisons directory: %s" %(output_dirname)
            misc_utils.make_dir(output_dirname)
	ht.output_samples_comparison(sample1_dirname,
                                     sample2_dirname,
                                     output_dirname,
                                     sample_labels=options.comparison_labels,
                                     use_compressed=use_compressed)
    ##
    ## Main interface based on SAM files
    ##
    if options.compute_genes_from_file != None:
        # Run on events given by file
        run_compute_genes_from_file(options)
    if options.compute_gene_psi != None:
        run_compute_gene_psi(options)
        
    ##
    ## Summarizing samples
    ##
    if options.summarize_samples:
	samples_dir = \
            os.path.abspath(os.path.expanduser(options.summarize_samples[0]))
        if options.summary_label != None:
            samples_label = options.summary_label
            print "Using summary label: %s" %(samples_label)
        else:
            samples_label = \
                os.path.basename(os.path.expanduser(samples_dir))
	assert(len(samples_label) >= 1)
	summary_output_dir = \
            os.path.abspath(os.path.join(os.path.expanduser(options.summarize_samples[1]),
                                         'summary'))
	if not os.path.isdir(summary_output_dir):
	    os.makedirs(summary_output_dir)
	    
	summary_filename = os.path.join(summary_output_dir,
					'%s.miso_summary' %(samples_label))
	summarize_sampler_results(samples_dir, summary_filename,
                                  use_compressed=use_compressed)

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

if __name__ == '__main__':
    main()
