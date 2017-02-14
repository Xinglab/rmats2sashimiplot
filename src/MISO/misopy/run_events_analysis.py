##
## Events analysis
##
import os
import csv
import sys
import subprocess
from collections import defaultdict

import pysam

import misopy
import misopy.gff_utils as gff_utils
import misopy.as_events as as_events
import misopy.run_miso as run_miso
import misopy.misc_utils as misc_utils
import misopy.exon_utils as exon_utils
from misopy.parse_csv import *
from misopy.settings import Settings, load_settings
from misopy.settings import miso_path as miso_settings_path
import misopy.cluster_utils as cluster_utils

miso_path = os.path.dirname(os.path.abspath(__file__))
manual_url = "http://genes.mit.edu/burgelab/miso/docs/"


def get_ids_passing_filter(gff_index_dir,
                           bam_filename,
                           output_dir):
    """
    Apply filter to events using bedtools and return
    only the events that meet the filter.
    """
    min_reads = 20
    settings = Settings.get()
    min_event_reads = Settings.get_min_event_reads()
    
    # Check that this was indexed with a version that outputs
    # genes.gff file
    genes_gff_fname = os.path.join(gff_index_dir,
                                   "genes.gff")
    if not os.path.isfile(genes_gff_fname):
        print "WARNING: Could not find \'genes.gff\' in %s - " \
              "skipping prefilter stage. Please reindex your " \
              "GFF file with the latest version to enable " \
              "prefiltering." %(gff_index_dir)
        return None
    print "Prefiltering reads..."
    coverage_fname = exon_utils.get_bam_gff_coverage(bam_filename,
                                                     genes_gff_fname,
                                                     output_dir)
    ids_passing_filter = []
    with open(coverage_fname) as coverage_in:
        for line in coverage_in:
            # Skip comments
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            # Get the counts field and the event ID
            # if it passes the filter
            counts = int(fields[9])
            if counts < min_event_reads:
                continue
            attribs = gff_utils.parse_gff_attribs(fields[8])
            if "ID" not in attribs:
                print "WARNING: No ID= found for line:\n%s\nSkipping..." \
                    %(line)
                continue
            event_id = attribs["ID"]
            ids_passing_filter.append(event_id)
    return ids_passing_filter
            

def check_gff_and_bam(gff_dir, bam_filename, main_logger,
                      num_genes=10000,
                      num_reads=10000,
                      given_read_len=None):
    """
    Look for chromosome headers mismatches between input GFF
    annotation and BAM filename.  Warn users if there are
    headers mismatches.
    """
    print "Checking your GFF annotation and BAM for mismatches..."
    # Check that BAM exists
    if not os.path.isfile(bam_filename):
        print "Error: BAM %s cannot be found." %(bam_filename)
        return
    # Check that a BAM header is available
    bam_index_fname = "%s.bai" %(bam_filename)
    if not os.path.isfile(bam_index_fname):
        main_logger.warning("Expected BAM index file %s not found." \
                            %(bam_index_fname))
        main_logger.warning("Are you sure your BAM file is indexed?")
    print "Checking if BAM has mixed read lengths..."
    bam_file = pysam.Samfile(bam_filename, "rb")
    n = 0
    seq_lens = {}
    for bam_read in bam_file:
        # Check more of the reads for mixed read lengths
        if n >= (num_reads * 10):
            break
        seq_lens[len(bam_read.seq)] = True
        n += 1
    all_seq_lens = seq_lens.keys()
    if len(all_seq_lens) > 1:
        msg = "Found mixed length reads in your BAM file: %s\n" \
              "MISO does not support mixed read lengths. Please " \
              "trim your reads to a uniform length and re-run " \
              "or run MISO separately on each read length. " \
              "Proceeding anyway, though this may " \
              "result in errors or inability to match reads to " \
              "your annotation.\n" %(bam_filename)
        msg += "Read lengths were: %s\n" %(",".join(map(str, all_seq_lens)))
        main_logger.warning(msg)
        time.sleep(5)
    else:
        print "Found reads of length %d in BAM." %(all_seq_lens[0])
        # Check the BAM read length against the read length that was
        # given
        if given_read_len != None:
            if all_seq_lens[0] != given_read_len:
                e = "Error: The given read length (%d) does not match " \
                    "the read length found in BAM (%d). Please re-run " \
                    "and pass the correct --read-len argument. " \
                    "Note that MISO does not support mixed read lengths " \
                    "in the same BAM file." \
                    %(given_read_len,
                      all_seq_lens[0])
                main_logger.error(e)
                sys.exit(1)
        time.sleep(5)
        
    genes_fname = os.path.join(gff_dir, "genes.gff")
    if not os.path.isfile(genes_fname):
        # No genes.gff found - warn user and abort headers check
        main_logger.warning("No genes.gff file found in %s. Did you index " \
                            "your GFF with an older version of MISO?" \
                            %(gff_dir))
        return
    # Read first few genes chromosomes
    gff_chroms = {}
    gff_starts_with_chr = False
    with open(genes_fname) as genes:
        n = 0
        for line in genes:
            if n >= num_genes:
                break
            curr_gff_chrom = line.strip().split("\t")[0]
            gff_chroms[curr_gff_chrom] = True
            # Record that GFF chroms start with chr if they do
            if curr_gff_chrom.startswith("chr"):
                gff_starts_with_chr = True
            n += 1
    # Read first few BAM reads chromosomes
    bam_file = pysam.Samfile(bam_filename, "rb")
    bam_chroms = {}
    bam_starts_with_chr = False
    n = 0
    for bam_read in bam_file:
        if n >= num_reads:
            break
        read_chrom = bam_file.getrname(bam_read.tid)
        bam_chroms[read_chrom] = True
        if str(read_chrom).startswith("chr"):
            bam_starts_with_chr = True
        n += 1
    mismatch_found = False
    if bam_starts_with_chr != gff_starts_with_chr:
        mismatch_found = True
    if mismatch_found:
        main_logger.warning("It looks like your GFF annotation file and your BAM " \
                            "file might not have matching headers (chromosome names.) " \
                            "If this is the case, your run will fail as no reads from " \
                            "the BAM could be matched up with your annotation.")
        main_logger.warning("Please see:\n\t%s\n for more information." %(manual_url))
        # Default: assume BAM starts with chr headers
        chr_containing = "BAM file (%s)" %(bam_filename)
        not_chr_containing = "GFF annotation (%s)" %(gff_dir)
        if not bam_starts_with_chr:
            # BAM does not start with chr, GFF does
            chr_containing, not_chr_containing = \
                not_chr_containing, chr_containing
        main_logger.warning("It looks like your %s contains \'chr\' chromosomes (UCSC-style) " \
                            "while your %s does not." %(chr_containing,
                                                        not_chr_containing))
        main_logger.warning("The first few BAM chromosomes were: %s" \
                            %(",".join(bam_chroms.keys())))
        print "BAM references: "
        print bam_file.references
        main_logger.warning("The first few GFF chromosomes were: %s" \
                            %(",".join(gff_chroms.keys())))
        main_logger.warning("Run is likely to fail or produce empty output. Proceeding " \
                            "anyway...")
        time.sleep(15)

        
def compute_psi(sample_filenames, output_dir, event_type,
                read_len, overhang_len,
                use_cluster=False,
                chunk_jobs=False,
                filter_events=True,
                events_info_filename=None,
                settings_filename=None):
    """
    Compute Psi values for skipped exons.  Sample filenames is a mapping from
    sample label to sample.

      - sample_filenames = [[sample_label1, sample_filename1],
                            [sample_label2, sample_filename2]]
      - output_dir: output directory
      - event_type: 'SE', 'RI', etc.
    """
    misc_utils.make_dir(output_dir)
    
    output_dir = os.path.join(output_dir, event_type)
    output_dir = os.path.abspath(output_dir)

    misc_utils.make_dir(output_dir)
	
    print "Computing Psi for events of type %s" %(event_type)
    print "  - samples used: ", sample_filenames.keys()

    for sample_label, sample_filename in sample_filenames.iteritems():
	print "Processing sample: label=%s, filename=%s" \
            %(sample_label, sample_filename)
	results_output_dir = os.path.join(output_dir, sample_label)
        misc_utils.make_dir(results_output_dir)

	# Load the set of counts and serialize them into JSON
	events = \
            as_events.load_event_counts(sample_filename,
                                        event_type,
                                        events_info_filename=events_info_filename)

	# Filter events
	if filter_events:
	    print "Filtering events..."
	    events.filter_events(settings=Settings.get())

	print "Running on a total of %d events." %(len(events.events))
	    
	events_filename = events.output_file(results_output_dir,
                                             sample_label)
	
	# Run MISO on them
	miso_cmd = "python %s --compute-two-iso-psi %s %s --event-type %s " \
                   "--read-len %d --overhang-len %d " \
                   %(os.path.join(miso_path, 'run_miso.py'),
                     events_filename,
                     results_output_dir,
                     event_type,
                     read_len,
                     overhang_len)
	if use_cluster:
	    if chunk_jobs:
		miso_cmd += ' --use-cluster --chunk-jobs %d' %(chunk_jobs)
	    else:
		miso_cmd += ' --use-cluster'
        print "Executing: %s" %(miso_cmd)
	if use_cluster:
	    print " - Using cluster"
	os.system(miso_cmd)


def greeting(parser=None):
    print "MISO (Mixture of Isoforms model)"
    print "Probabilistic analysis of RNA-Seq data to detect " \
          "differential isoforms"
    print "Use --help argument to view options.\n"
    if parser is not None:
        parser.print_help()


def main():
    print "MISO (Mixture of Isoforms model)"
    print "To run MISO, please use \"miso\" instead."
		    
if __name__ == '__main__':
    main()
    
