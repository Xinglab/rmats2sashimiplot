# -*- mode: python; -*-
##
## Utilities for working with paired-end reads and
## fragment distributions
##

import os
import sys
import glob
import time
import pysam

from scipy import *
from numpy import *

import misopy
from misopy.parse_csv import *
import misopy.Gene as gene_utils
import misopy.sam_utils as sam_utils
import misopy.exon_utils as exon_utils
from collections import defaultdict

def get_insert_dist_array(interval_to_paired_dists,
                          delim='\t'):
    """
    Read insert length distribution as array of numbers.
    """
    insert_dist = []
    for interval, paired_dists in interval_to_paired_dists.iteritems():
        insert_dist.extend(paired_dists)
    return array(insert_dist)


def parse_insert_len_params(insert_len_header):
    params = {}
    insert_len_header = insert_len_header.strip()
    if insert_len_header.startswith("#"):
        insert_len_header = insert_len_header[1:]
    for param in insert_len_header.split(","):
        p, v = param.split("=")
        params[p] = v
    return params

        
def filter_insert_len(interval_to_dists,
                     sd_max):
    # Get vector of insert lengths 
    insert_dist = get_insert_dist_array(interval_to_dists)
    
    mu, sdev, dispersion, num_pairs = \
        compute_insert_len_stats(insert_dist)

    filtered_interval_to_dists = defaultdict(list)
    
    min_cutoff = mu - (sd_max * sdev)
    max_cutoff = mu + (sd_max * sdev)
    print "Excluding values < %.2f or > %.2f" \
          %(min_cutoff, max_cutoff)

    for interval, dists in interval_to_dists.iteritems():
        dists = array(dists)
        filtered_dists = delete(dists, nonzero(dists < min_cutoff)[0])
        filtered_dists = delete(filtered_dists,
                                nonzero(dists > max_cutoff)[0])
        filtered_interval_to_dists[interval] = filtered_dists
                                             
    return filtered_interval_to_dists
               

def load_insert_len(insert_dist_filename,
                    delim='\t'):
    print "Loading insert length from: %s" %(insert_dist_filename)
    insert_dist_file = open(insert_dist_filename, "r")
    insert_lens = []
    params_header = insert_dist_file.readline().strip()
    
    # Get parameters of distribution from header
    params = parse_insert_len_params(params_header)
    
    for line in insert_dist_file:
        # Skip header
        if line.startswith("#"):
            continue
        insert_fields = line.strip().split(delim)
        if len(insert_fields) != 2:
            continue
        lens_list = insert_fields[1].split(",")
        curr_lens = [int(l) for l in lens_list]
        insert_lens.extend(curr_lens)
    insert_dist = array(insert_lens)
    insert_dist_file.close()
    return insert_dist, params


def bedtools_map_bam_to_bed(bam_filename, gff_intervals_filename):
    """
    Map BAM file to GFF intervals and return the result as a
    BED file.

    Returns a stream to a BED file with the results
    """
    bedtools_cmd = "intersectBed -abam %s -b %s -wa -wb -bed -f 1" \
                   %(bam_filename, gff_intervals_filename)

    print "Executing: %s" %(bedtools_cmd)

    if (not os.path.isfile(bam_filename)) or \
       (not os.path.isfile(gff_intervals_filename)):
        raise Exception, "Error: %s or %s do not exist." \
            %(bam_filename,
              gff_intervals_filename)
    bed_stream = os.popen(bedtools_cmd)
    return bed_stream


def parse_tagBam_intervals(bam_read,
                           gff_coords=True):
    """
    Return a list of intervals that are present in the current
    BAM line returned by tagBam. These intervals are encoded
    in the YB option of the BAM flag.

    - If convert_coords is True, we add 1 to the
      BAM coordinate to make it 1-based
    """
    gff_aligned_regions = bam_read.opt("YB")
    parsed_regions = gff_aligned_regions.split("gff:")[1:]

    gff_intervals = []

    for region in parsed_regions:
        strand = region.split(",")[3]
        chrom, coord_field = region.split(",")[0].split(":")
        region_start, region_end = coord_field.split("-")
        region_start, region_end = int(region_start), \
                                   int(region_end)
        if gff_coords:
            region_start += 1
        curr_interval_str = "%s:%d-%d:%s" \
                            %(chrom,
                              region_start,
                              region_end,
                              strand)
        gff_intervals.append(curr_interval_str)
    return gff_intervals


def compute_inserts_from_paired_mates(paired_reads):
    """
    Get insert lengths from paired-up paired ends reads
    aligned to a set of constitutive exon intervals.

    Return mapping from intervals to distances of read pairs
    that land in them.
    """
    # Mapping from interval to 
    interval_to_paired_dists = defaultdict(list)
    num_skipped = 0
    num_kept = 0
    for read_id, read_pair in paired_reads.iteritems():
        to_skip = False
        # Get the intervals that each read pair lands in
        # Consider here only the mate pairs that map to
        # the same interval, and to exactly one interval, and
        # not in a junction
        left_mate, right_mate = read_pair
        left_mate_intervals = parse_tagBam_intervals(left_mate)
        right_mate_intervals = parse_tagBam_intervals(right_mate)
        
        # If either of the mates lands in more than one set of intervals,
        # discard it.
        if (len(left_mate_intervals) != 1 or \
            len(right_mate_intervals) != 1):
            to_skip = True
        elif left_mate_intervals[0] != right_mate_intervals[0]:
            # If each maps to one interval, but it's not the same,
            # also discard it.
            to_skip = True
        elif (len(left_mate.cigar) != 1 or \
              len(right_mate.cigar) != 1):
            # One of the read mates was in a junction
            to_skip = True
        elif (left_mate.cigar[0][0] != 0 or \
              right_mate.cigar[0][0] != 0):
            # Both CIGAR operations must be M (matches)
            to_skip = True

        if to_skip:
            # One of the conditions was violated, so skip read pair
            num_skipped += 1
            continue

        # We have a match, so compute insert length distance,
        # defined as the distance between the start position
        # of the left and the end position of the right mate
        left_start = left_mate.pos
        left_end = sam_utils.cigar_to_end_coord(left_start,
                                                left_mate.cigar)

        right_start = right_mate.pos
        right_end = sam_utils.cigar_to_end_coord(right_start,
                                                 right_mate.cigar)

        # Get the current GFF interval string
        curr_gff_interval = left_mate_intervals[0]
        
        # Insert length is right.end - left.start + 1
        insert_len = right_end - left_start + 1

        if insert_len <= 0:
            print "WARNING: 0 or negative insert length detected " \
                  "in region %s." %(curr_gff_interval)
            continue
        
        interval_to_paired_dists[curr_gff_interval].append(insert_len)
        num_kept += 1

    print "Used %d paired mates, threw out %d" \
          %(num_kept, num_skipped)

    return interval_to_paired_dists
            
    
def compute_insert_len(bams_to_process,
                       const_exons_gff_filename,
                       output_dir,
                       min_exon_size,
                       no_bam_filter=False,
                       sd_max=2):
    """
    Compute insert length distribution and output it to the given
    directory.

    Arguments:

    - bams_to_process: a list of BAM files to process
    - const_gff_filename: GFF with constitutive exons
    """
    bams_str = "\n  ".join(bams_to_process)
    num_bams = len(bams_to_process)
    print "Computing insert length distribution of %d files:\n  %s" \
          %(num_bams, bams_str)
    print "  - Using const. exons from: %s" %(const_exons_gff_filename)
    print "  - Outputting to: %s" %(output_dir)
    print "  - Minimum exon size used: %d" %(min_exon_size)

    if not os.path.isdir(output_dir):
        print "Making directory: %s" %(output_dir)
        os.makedirs(output_dir)

    all_constitutive = True

    const_exons, f = \
        exon_utils.get_const_exons_by_gene(const_exons_gff_filename,
                                           output_dir,
                                           # Treat all exons as constitutive
                                           all_constitutive=True,
                                           min_size=min_exon_size)
    filter_reads = not no_bam_filter

    if filter_reads:
        print "Filtering BAM reads"
    else:
        print "Turning off filtering of BAM reads"
        
    for bam_filename in bams_to_process:
        t1 = time.time()
        output_filename = os.path.join(output_dir,
                                       "%s.insert_len" \
                                       %(os.path.basename(bam_filename)))
        if not os.path.isfile(bam_filename):
            print "Cannot find BAM file %s" %(bam_filename)
            print "Quitting..."
            sys.exit(1)
        print "Fetching reads in constitutive exons"
        mapped_bam_filename = exon_utils.map_bam2gff(bam_filename,
                                                     const_exons_gff_filename,
                                                     output_dir)
        if mapped_bam_filename == None:
            raise Exception, "Error: Insert length computation failed."

        # Load mapped BAM filename
        mapped_bam = pysam.Samfile(mapped_bam_filename, "rb")
        ###
        ### TODO: Rewrite this so that you only pair reads within an interval
        ###
        paired_reads = sam_utils.pair_sam_reads(mapped_bam,
                                                filter_reads=filter_reads)
        num_paired_reads = len(paired_reads)

        if num_paired_reads == 0:
            print "WARNING: no paired mates in %s. Skipping...\n"\
                  "Are you sure the read IDs match? If your BAM paired flags are "\
                  "unset, try using --no-bam-filter." \
                  %(bam_filename)
            continue
        print "Using %d paired mates" %(num_paired_reads)
        interval_to_paired_dists = compute_inserts_from_paired_mates(paired_reads)
        summarize_insert_len_dist(interval_to_paired_dists, output_filename,
                                  sd_max=sd_max)
        t2 = time.time()
        print "Insert length computation took %.2f seconds." %(t2 - t1)


# def pair_reads_from_bed_intervals(bed_stream):
#     """
#     Match up read mates with each other, indexed by the BED interval
#     that they fall in.

#     Return a dictionary of BED region mapping to a set of read pairs.

#     Arguments:

#     - bed_filename: file with BED reads and the region they map to.

#     Returns. 
#     """
#     return

# def compute_insert_len(bam_filename, gff_filename, output_dir,
#                        min_exon_size):
#     """
#     Compute insert length distribution and output it to the given
#     directory.
#     """
#     print "Computing insert length distribution of %s" %(bam_filename)
#     print "  - Using gene models from: %s" %(gff_filename)
#     print "  - Outputting to: %s" %(output_dir)
#     print "  - Minimum exon size used: %d" %(min_exon_size)

#     if not os.path.isdir(output_dir):
#         print "Making directory: %s" %(output_dir)
#         os.makedirs(output_dir)

#     output_filename = os.path.join(output_dir,
#                                    "%s.insert_len" %(os.path.basename(bam_filename)))

#     # Load BAM file with reads
#     bamfile = sam_utils.load_bam_reads(bam_filename)
    
#     # Load the genes from the GFF
#     print "Loading genes from GFF..."
#     t1 = time.time()
#     gff_genes = gene_utils.load_genes_from_gff(gff_filename)
#     t2 = time.time()
#     print "  - Loading genes from GFF took %.2f seconds" %(t2 - t1)

#     insert_lengths = []

#     t1 = time.time()

#     relevant_region = 0
    
#     for gene_id, gene_info in gff_genes.iteritems():
#         gene_obj = gene_info["gene_object"]

#         # Get all the constitutive parts
#         const_parts = gene_obj.get_const_parts()

#         chrom = gene_obj.chrom

#         # Consider only the large constitutive parts
#         for part in const_parts:
#             if part.len >= min_exon_size:
#                 # Get all the reads that land in the coordinates of the exon
#                 try:
#                     exon_reads = bamfile.fetch(chrom, part.start, part.end)
#                 except ValueError:
#                     print "Could not fetch from region: ", chrom, part.start, part.end
#                     continue

#                 # Pair all the paired-end reads that land there
#                 paired_reads = sam_utils.pair_sam_reads(exon_reads)
#                 num_paired_reads = len(paired_reads)

#                 if num_paired_reads == 0:
#                     continue

#                 print "Found %d region" %(relevant_region)
#                 relevant_region += 1

#                 # Compute the insert length of each read
#                 for read_pair_id, read_pair in paired_reads.iteritems():
#                     if len(read_pair) != 2:
#                         # Skip non-paired reads
#                         continue
                    
#                     left_read, right_read = read_pair
#                     insert_len = right_read.pos - left_read.pos + 1

#                     if insert_len > 0:
#                         insert_lengths.append(insert_len)
#                     else:
#                         print "Negative or zero insert length ignored..."

#     # Output results to file
#     output_file = open(output_filename, 'w')
#     insert_length_str = "\n".join(map(str, insert_lengths))
#     output_file.write(insert_length_str)
#     output_file.close()
                    
#     t2 = time.time()
#     print "Insert length computation took %.2f seconds." %(t2 - t1)


def output_insert_len_dist(interval_to_paired_dists,
                           output_file):
    """
    Output insert length distribution indexed by regions.
    """
    header = "#%s\t%s\n" %("region", "insert_len")
    output_file.write(header)

    for region, insert_lens in interval_to_paired_dists.iteritems():
        if len(insert_lens) == 0:
            continue
        str_lens = ",".join([str(l) for l in insert_lens])
        output_line = "%s\t%s\n" %(region, str_lens)
        output_file.write(output_line)


def compute_insert_len_stats(insert_dist):
    """
    Return insert length statistics.
    """
    # Compute mean and standard deviation of insert
    # length distribution
    mu = mean(insert_dist)
    sdev = std(insert_dist)
    
    # Compute dispersion (d), where
    #
    # d = sdev / sqrt(mean)
    #
    # dispersion measures how variable
    # the insert length distribution is
    # about the mean
    dispersion = sdev / sqrt(float(mu))
    
    # Number of read pairs used 
    num_pairs = len(insert_dist)

    return mu, sdev, dispersion, num_pairs
    

def summarize_insert_len_dist(interval_to_paired_dists,
                              output_filename,
                              sd_max=2):
    """
    Summarize insert len distributions.
    """
    print "Summarizing insert length distribution.."
    print "  - Output file: %s" %(output_filename)

    output_file = open(output_filename, "w")

    print "Removing values %d-many deviations outside the mean" \
          %(sd_max)

    # Filter insert length distribution based on sd_max
    filtered_interval_to_dist = filter_insert_len(interval_to_paired_dists,
                                                  sd_max)
    filtered_insert_dist = get_insert_dist_array(filtered_interval_to_dist)

    if len(filtered_insert_dist) == 0:
        print "Error: Could not find any properly mated pairs to " \
              "compute insert length with. Are you sure your BAM reads " \
              "are properly paired and map the chromosome headers in the " \
              "constitutive exon file?"
        sys.exit(1)

    mu, sdev, dispersion, num_pairs = \
        compute_insert_len_stats(filtered_insert_dist)
        
    print "mean\tsdev\tdispersion"
    print "%.1f\t%.1f\t%.1f" \
          %(mu, sdev, dispersion)
    min_insert = min(filtered_insert_dist)
    max_insert = max(filtered_insert_dist)
    print "min insert: %d" %(min_insert)
    print "max insert: %d" %(max_insert)

    # Write headers
    header_line = "#%s=%.1f,%s=%.1f,%s=%.1f,%s=%d\n" \
                  %("mean", mu,
                    "sdev", sdev,
                    "dispersion", dispersion,
                    "num_pairs", num_pairs)
    output_file.write(header_line)

    # Write raw insert lengths indexed by region
    output_insert_len_dist(filtered_interval_to_dist,
                           output_file)
    output_file.close()


def greeting():
    print "Utility for computing insert length distributions from paired-end " \
          "BAM files."
    print "Part of MISO (Mixture of Isoforms model)\n"    
    print "See --help for usage.\n"


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--compute-insert-len", dest="compute_insert_len", nargs=2, default=None,
                      help="Compute insert length for given sample. Takes as input "
                      "(1) a comma-separated list of sorted, indexed BAM files with headers "
                      "(or a single BAM filename), (2) a GFF file with constitutive exons. "
                      "Outputs the insert length distribution into the output directory.")
    parser.add_option("--no-bam-filter", dest="no_bam_filter", action="store_true", default=False,
                      help="If provided, this ignores the BAM file flags that state whether the read was paired "
                      "or not, and instead uses only the read IDs to pair up the mates. Use this if your "
                      "paired-end BAM was the result of a samtools merge operation.")
    parser.add_option("--min-exon-size", dest="min_exon_size", nargs=1, type="int", default=500,
                      help="Minimum size of constitutive exon (in nucleotides) that should be used "
                      "in the computation. Default is 500 bp.")
    parser.add_option("--sd-max", dest="sd_max", nargs=1, default=2, type="int",
                      help="Number of standard deviations used to define outliers. By default, set "
                      "to 2, meaning that any points at least 2*sigma away from the mean of the "
                      "insert length distribution will be discarded.")
    parser.add_option("--output-dir", dest="output_dir", nargs=1, default=None,
                      help="Output directory.")
    (options, args) = parser.parse_args()

    if options.compute_insert_len is None:
        greeting()
        return

    if options.output_dir is None:
        greeting()
        
        print "Error: need --output-dir."
        return
        
    output_dir = os.path.abspath(os.path.expanduser(options.output_dir))

    sd_max = options.sd_max

    if options.compute_insert_len != None:
        bams_to_process = [os.path.abspath(os.path.expanduser(f)) for f in \
                           options.compute_insert_len[0].split(",")]
        gff_filename = os.path.abspath(os.path.expanduser(options.compute_insert_len[1]))
        compute_insert_len(bams_to_process, gff_filename, output_dir,
                           options.min_exon_size,
                           no_bam_filter=options.no_bam_filter,
                           sd_max=sd_max)

if __name__ == "__main__":
    main()
