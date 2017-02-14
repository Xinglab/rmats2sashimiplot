##
## Utilities for handling SAM/BAM reads
##

from collections import defaultdict

import misopy
from misopy.Gene import load_genes_from_gff

import os
import time
import pysam
import binascii
import ctypes

from numpy import array
from scipy import *

def cigar_to_end_coord(start, cigar):
    """
    Compute the end coordinate based on the CIGAR string.

    Assumes the coordinate is 1-based.
    """
    #print start, cigar
    # Parse cigar part
    #for cigar_part in cigar:
    #    cigar_type, seq_len = cigar_part
    #    offset += seq_len
    offset = sum([cigar_part[1] for cigar_part in cigar])
    end = start + offset - 1
    return end

def single_end_read_to_isoforms(read, gene, read_len, overhang_len=1):
    """
    Align single-end SAM read to gene's isoforms.
    """
    start = read.pos + 1
    end = cigar_to_end_coord(start, read.cigar)    
    
    assert(start < end)

    alignment, isoform_coords = gene.align_read_to_isoforms_with_cigar(read.cigar, start, end, read_len,
                                                                       overhang_len)
    return alignment


def paired_read_to_isoforms(paired_read, gene, read_len,
                            overhang_len=1):
    """
    Align paired-end SAM read to gene's isoforms.
    """
    left_read = paired_read[0]
    right_read = paired_read[1]

    # Convert to 1-based coordinates
    left_start = left_read.pos + 1
    right_start = right_read.pos + 1

    # Get end coordinates of each read
    left_end = cigar_to_end_coord(left_start, left_read.cigar)

    assert(left_start < left_end)
    
    right_end = cigar_to_end_coord(right_start, right_read.cigar)

    assert(right_start < right_end)

    alignment = None
    frag_lens = None
    
    # Throw out reads that posit a zero or less than zero fragment
    # length
#    print "LEFT read start,end: ", left_start, left_end
#    print "RIGHT read start,end: ", right_start, right_end
    if left_start > right_start:
        return None, None
    else:
        alignment, frag_lens = gene.align_read_pair_with_cigar(
            left_read.cigar, left_start, left_end, 
            right_read.cigar, right_start, right_end,
            read_len=read_len, overhang=overhang_len)
        
    
#    assert(left_start < right_start), "Reads scrambled?"
#    assert(left_end < right_start), "Reads scrambled: left=(%d, %d), right=(%d, %d)" \
#                    %(left_start, left_end, right_start, right_end)
    
#    print "Read: ", (left_start, left_end), " - ", (right_start, right_end)
#    print "  - Alignment: ", alignment, " frag_lens: ", frag_lens
#    print "  - Sequences: "
#    print "  - %s\t%s" %(left_read.seq, right_read.seq)
    
    return alignment, frag_lens

# def paired_read_to_isoforms(paired_read, gene, read_len=36,
#                             overhang_len=1):
#     """
#     Align paired-end SAM read to gene's isoforms.
#     """
#     left_read = paired_read[0]
#     right_read = paired_read[1]

#     # Convert to 1-based coordinates
#     left_start = left_read.pos + 1
#     right_start = right_read.pos + 1

#     # Get end coordinates of each read
#     left_end = cigar_to_end_coord(left_start, left_read.cigar)

#     assert(left_start < left_end)
    
#     right_end = cigar_to_end_coord(right_start, right_read.cigar)

#     assert(right_start < right_end)

#     alignment = None
#     frag_lens = None
    
#     # Throw out reads that posit a zero or less than zero fragment
#     # length
#     if left_start > right_start or left_end > right_start:
#         return None, None
#     else:
#         alignment, frag_lens = gene.align_read_pair(left_start, left_end,
#                                                     right_start, right_end,
#                                                     read_len=read_len,
#                                                     overhang=overhang_len)
        
    
# #    assert(left_start < right_start), "Reads scrambled?"
# #    assert(left_end < right_start), "Reads scrambled: left=(%d, %d), right=(%d, %d)" \
# #                    %(left_start, left_end, right_start, right_end)
    
# #    print "Read: ", (left_start, left_end), " - ", (right_start, right_end)
# #    print "  - Alignment: ", alignment, " frag_lens: ", frag_lens
# #    print "  - Sequences: "
# #    print "  - %s\t%s" %(left_read.seq, right_read.seq)
    
#     return alignment, frag_lens
                         
def load_bam_reads(bam_filename,
                   template=None):
    """
    Load a set of indexed BAM reads.
    """
    print "Loading BAM filename from: %s" %(bam_filename)
    bam_filename = os.path.abspath(os.path.expanduser(bam_filename))
    bamfile = pysam.Samfile(bam_filename, "rb",
                            template=template)
    return bamfile


def fetch_bam_reads_in_gene(bamfile, chrom, start, end,
                            gene=None):
    """
    Align BAM reads to the gene model.
    """
    gene_reads = []

    if chrom in bamfile.references:
        pass
    else:
        chrom_parts = chrom.split("chr")
        if len(chrom_parts) <= 1:
            chrom = chrom_parts[0]
        else:
            chrom = chrom_parts[1]

    try:
        gene_reads = bamfile.fetch(chrom, start, end)
    except ValueError:
        print "Cannot fetch reads in region: %s:%d-%d" %(chrom,
                                                         start,
                                                         end)
    except AssertionError:
        print "AssertionError in region: %s:%d-%d" %(chrom,
                                                     start,
                                                     end)
        print "  - Check that your BAM file is indexed!"
    return gene_reads


def flag_to_strand(flag):
    """
    Takes integer flag as argument.
    Returns strand ('+' or '-') from flag.
    """
    if flag & 16:
        return "-"
    return "+"


def strip_mate_id(read_name):
    """
    Strip canonical mate IDs for paired end reads, e.g.
    
    #1, #2

    or:

    /1, /2
    """
    if read_name.endswith("/1") or read_name.endswith("/2") or \
       read_name.endswith("#1") or read_name.endswith("#2"):
        read_name = read_name[0:-3]
    return read_name

    
def pair_sam_reads(samfile,
                   filter_reads=True,
                   return_unpaired=False,
                   strand_rule=None):
    """
    Pair reads from a SAM file together.
    """
    paired_reads = defaultdict(list)
    unpaired_reads = {}

    for read in samfile:
        curr_name = read.qname

        # Strip canonical mate IDs 
        curr_name = strip_mate_id(curr_name)
        
        if filter_reads:
            # Skip reads that failed QC or are unmapped
            if read.is_qcfail or read.is_unmapped or \
               read.mate_is_unmapped or (not read.is_paired):
                unpaired_reads[curr_name] = read
                continue
        paired_reads[curr_name].append(read)
        # Ensure that the reads that were paired are
        # in the right order - i.e., that read1 is
        # first and read2 follows.
        if len(paired_reads[curr_name]) == 2:
            if strand_rule == "fr-firststrand":
                # Thanks to Renee Sears:
                # For fr-firststrand the /1 read should only be left of the
                # right if it is on the '+' strand whereas the /2 read
                # should be to the left if it is on the '+' strand
                if paired_reads[curr_name][0].is_read1 and \
                   paired_reads[curr_name][0].is_reverse:
                    paired_reads[curr_name] = paired_reads[curr_name][::-1]
                if paired_reads[curr_name][0].is_read2 and \
                    paired_reads[curr_name][0].is_reverse:
                    paired_reads[curr_name] = paired_reads[curr_name][::-1]
                    
    to_delete = []
    num_unpaired = 0
    num_total = 0

    for read_name, read in paired_reads.iteritems():
        if len(read) != 2:
            unpaired_reads[read_name] = read
            num_unpaired += 1
            # Delete unpaired reads
            to_delete.append(read_name)
            continue
        left_read, right_read = read[0], read[1]

        # Check that read mates are on opposite strands
        left_strand = flag_to_strand(left_read.flag)
        right_strand = flag_to_strand(right_read.flag)

        if left_strand == right_strand:
            # Skip read pairs that are on the same strand
            to_delete.append(read_name)
            continue
        
        if left_read.pos > right_read.pos:
            print "WARNING: %s left mate starts later than right "\
                  "mate" %(left_read.qname)
        num_total += 1

    # Delete reads that are on the same strand
    for del_key in to_delete:
        del paired_reads[del_key]

    print "Filtered out %d read pairs that were on same strand." \
        %(len(to_delete))
    print "Filtered out %d reads that had no paired mate." \
        %(num_unpaired)
    print "  - Total read pairs: %d" %(num_total)

    if not return_unpaired:
        return paired_reads
    else:
        return paired_reads, unpaired_reads


# Global variable containing CIGAR types for conversion
CIGAR_TYPES = ('M', 'I', 'D', 'N', 'S', 'H', 'P')

def sam_cigar_to_str(sam_cigar):
    """
    Convert pysam CIGAR list to string format.
    """
    # First element in sam CIGAR list is the CIGAR type
    # (e.g. match or insertion) and the second is
    # the number of nucleotides
    #cigar_str = "".join(["%d%s" %(c[1], CIGAR_TYPES[c[0]]) \
    #                     for c in sam_cigar])
    #### OPTIMIZED VERSION
    cigar_str = ""
    if sam_cigar is None:
        return cigar_str
    for c in sam_cigar:
        cigar_str += "%d%s" %(c[1], CIGAR_TYPES[c[0]])
    return cigar_str


def read_matches_strand(read,
                        target_strand,
                        strand_rule,
                        paired_end=None):
    """
    Check if a read matches strand.

    - target_strand: the annotation strand ('+' or '-')
    - strand_rule: the strand rule, i.e.
      ('fr-unstranded', 'fr-firststrand', or 'fr-secondstrand')
    """
    if strand_rule == "fr-unstranded":
        return True
    if strand_rule == "fr-secondstrand":
        raise Exception, "fr-secondstrand currently unsupported."
    matches = False
    if paired_end is not None:
        ## Paired-end reads
        read1, read2 = read
        if strand_rule == "fr-firststrand":
            # fr-firststrand: means that the *first* of the mates
            # must match the strand. Superfluous in light of
            # switching in pair_sam_reads()
            if target_strand == "+":
                return (flag_to_strand(read1.flag) == "+")
            elif target_strand == "-":
                return (flag_to_strand(read2.flag) == "-")
        else:
            raise Exception, "Unknown strandedness rule."
    else:
        ## Single-end reads
        if strand_rule == "fr-firststrand":
            # fr-firststrand: We sequence the first read only, so it must
            # match the target strand
            matches = (flag_to_strand(read.flag) == target_strand)
        else:
            raise Exception, "Unknown strandedness rule."
    return matches 


def sam_parse_reads(samfile,
                    paired_end=False,
                    strand_rule=None,
                    target_strand=None,
                    given_read_len=None):
    """
    Parse the SAM reads. If paired-end, pair up the mates
    together.

    Also forces the strandedness convention, discarding
    reads that do not match the correct strand.

    Kwargs:
    - paired_end: whether paired-end or not
    - strand_rule: specifies the strandedness convention. Can be
      'fr-unstranded', 'fr-firststrand' or 'fr-secondstrand'.
    - target_strand: specifies the strand to match, i.e. the
      annotation strand. Can be '+' or '-'.
    - given_read_len: The read length given to MISO by the user.
    If passed to this function, it will filter out all reads 
    that do not have this length (e.g. in mixed read length
    BAM file.)
    """
    read_positions = []
    read_cigars = []
    num_reads = 0

    check_strand = True
    # Determine if we need to check strandedness of reads.
    # If we're given an unstranded convention, or if we're
    # not given a target strand, then assume that there's
    # no need to check strandedness.
    if (strand_rule is None) or \
       (strand_rule is "fr-unstranded") or \
       (target_strand is None):
        # No need to check strand
        check_strand = False

    # Track number of reads discarded due to strand
    # violations, if strand-specific
    num_strand_discarded = 0
    if paired_end:
        # Pair up the reads 
        paired_reads = pair_sam_reads(samfile,
                                      strand_rule=strand_rule)
        # Process reads into format required by fastmiso
        # MISO C engine requires pairs to follow each other in order.
        # Unpaired reads are not supported.
        for read_id, read_info in paired_reads.iteritems():
            if check_strand:
                # Check strand
                if not read_matches_strand(read_info,
                                           target_strand,
                                           strand_rule,
                                           paired_end=paired_end):
                    # Skip reads that don't match strand
                    num_strand_discarded += 1
                    continue
            read1, read2 = read_info
            if (read1.cigar is None) or (read2.cigar is None):
                continue
            # Filter on given read length here for PAIRED-END
            # If either mate is not of the given read length,
            # skip it
            if given_read_len is not None:
                if (read1.rlen != given_read_len) or \
                   (read2.rlen != given_read_len):
                   continue
            # Read positions and cigar strings are collected
            read_positions.append(int(read1.pos))
            read_positions.append(int(read2.pos))
            read_cigars.append(sam_cigar_to_str(read1.cigar))
            read_cigars.append(sam_cigar_to_str(read2.cigar))
            num_reads += 1
    else:
        # Single-end
        for read in samfile:
            if read.cigar is None:
                continue
            # Filter on given read length here for SINGLE-END
            if given_read_len is not None:
                # Skip reads that don't have the given read length
                if read.rlen != given_read_len:
                    continue
            if check_strand:
                if not read_matches_strand(read,
                                           target_strand,
                                           strand_rule,
                                           paired_end=paired_end):
                    # Skip reads that don't match strand
                    num_strand_discarded += 1
                    continue
            read_positions.append(int(read.pos))
            read_cigars.append(sam_cigar_to_str(read.cigar))
            num_reads += 1

    if check_strand:
        print "No. reads discarded due to strand violation: %d" \
            %(num_strand_discarded)

    reads = (tuple(read_positions),
             tuple(read_cigars))

    return reads, num_reads
    

def sam_pe_reads_to_isoforms(samfile, gene, read_len, overhang_len):
    """
    Align read pairs (from paired-end data set) to gene.

    Returns alignment of paired-end reads (with insert lengths)
    to gene and number of read pairs aligned.
    """
    paired_reads = pair_sam_reads(samfile)

    num_read_pairs = 0

    pe_reads = []

    k = 0

    for read_id, read_pair in paired_reads.iteritems():
        if len(read_pair) != 2:
            # Skip reads with no pair
            continue
        
        alignment, frag_lens = paired_read_to_isoforms(read_pair, gene,
                                                       read_len, overhang_len)
        
        # Skip reads that are not consistent with any isoform
        if any(array(alignment) == 1):
            pe_reads.append([alignment, frag_lens])
            num_read_pairs += 1
        else:
#            print "read %s inconsistent with all isoforms" %(read_id)
            k += 1

    print "Filtered out %d reads that were not consistent with any isoform" %(k)
    return pe_reads, num_read_pairs


def sam_se_reads_to_isoforms(samfile, gene, read_len,
                             overhang_len):
    """
    Align single-end reads to gene.
    """
    num_reads = 0

    alignments = []

    num_skipped = 0
    
    for read in samfile:
        alignment = single_end_read_to_isoforms(read, gene, read_len,
                                                overhang_len)
        if 1 in alignment:
            # If the read aligns to at least one of the isoforms, keep it
            alignments.append(alignment)
            num_reads += 1
        else:
            num_skipped += 1

    print "Skipped total of %d reads." %(num_skipped)
        
    return alignments, num_reads


def sam_reads_to_isoforms(samfile, gene, read_len, overhang_len,
                          paired_end=False):
    """
    Align BAM reads to the gene model.
    """
    print "Aligning reads to gene..."
    t1 = time.time()

    if paired_end != None:
        # Paired-end reads
        reads, num_reads = sam_pe_reads_to_isoforms(samfile, gene, read_len,
                                                    overhang_len)
    else:
        # Single-end reads
        reads, num_reads = sam_se_reads_to_isoforms(samfile, gene, read_len,
                                                    overhang_len)
    
    t2 = time.time()
    print "Alignment to gene took %.2f seconds (%d reads)." %((t2 - t1),
                                                              num_reads)
    return reads





def main():
    pass


if __name__ == "__main__":
    main()
