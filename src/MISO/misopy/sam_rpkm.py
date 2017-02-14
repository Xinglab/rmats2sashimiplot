#!/usr/bin/env python
##
## Compute RPKM from SAM/BAM files
##

import os
import time

import scipy
import numpy

from scipy import *
from numpy import *

import misopy
import misopy.sam_utils as sam_utils
from misopy.Gene import load_genes_from_gff
from misopy.parse_csv import *

import pysam


def rpkm_per_region(region_lens, region_counts, read_len,
                    num_total_reads):
    """
    Compute RPKM for the set of regions (defined by their lengths)
    and the counts in the region, assuming the given read length.
    """
    lens = array(region_lens)
    
    num_reads = sum(region_counts)

    # Number of mappable positions (in KB)
    num_positions_in_kb = (sum(lens - read_len + 1)) / 1e3

    # Number of reads (in millions)
    num_reads_in_millions = num_total_reads / 1e6

    # Reads counts 
    rpkm = (num_reads/num_positions_in_kb) / num_reads_in_millions

    return rpkm


class Counter:
    mCounts = 0
    def __call__(self, alignment):
        self.mCounts += 1

        
def count_total_reads(bam_filename):
    """
    Return total number of proper reads in BAM file.
    """
    bamfile = sam_utils.load_bam_reads(bam_filename)
    num_total_reads = 0

    for r in bamfile:
        num_total_reads += 1

    return num_total_reads


def compute_rpkm(gff_filename, bam_filename, read_len,
                 output_dir):
    """
    Compute RPKMs for genes listed in GFF based on BAM reads.
    """
    print "Computing RPKMs..."
    print "  - GFF filename: %s" %(gff_filename)
    print "  - BAM filename: %s" %(bam_filename)
    print "  - Output dir: %s" %(output_dir)
    
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    output_filename = os.path.join(output_dir,
                                   "%s.rpkm" %(os.path.basename(bam_filename)))
    print "Outputting RPKMs to: %s" %(output_filename)

    rpkm_fieldnames = ['gene_id', 'rpkm', 'const_exon_lens',
                       'num_reads']

    # Parse the GFF into genes
    print "Parsing GFF into genes..."
    t1 = time.time()
    gff_genes = load_genes_from_gff(gff_filename)
    t2 = time.time()
    print "Parsing took %.2f seconds" %(t2 - t1)

    # Load the BAM file
    bamfile = sam_utils.load_bam_reads(bam_filename)

    print "Counting all reads..."
    t1 = time.time()
    num_total_reads = count_total_reads(bam_filename)
    t2 = time.time()
    print "Took: %.2f seconds" %(t2 - t1)

    print "Number of total reads in BAM file: %d" %(num_total_reads)

    num_genes = 0

    rpkms_dictlist = []

    exons_too_small = {}
    num_no_const = 0

    for gene_id, gene_info in gff_genes.iteritems():
        # Get the gene object
        gene = gene_info['gene_object']

        # Get constitutive exons
        const_exons = gene.get_const_parts()

        num_reads = []
        exon_lens = []

        regions_counted = {}

        if not gene.chrom.startswith("chr"):
            chrom = "chr%s" %(gene.chrom)
        else:
            chrom = gene.chrom

        if "random" in chrom:
            print "Skipping random chromosome gene: %s, %s" \
                  %(gene_id, chrom)
            continue

        if len(const_exons) == 0:
            print "Gene %s has no constitutive regions!" %(gene_id)
            num_no_const += 1
            continue

        total_num_reads = 0

        for exon in const_exons:
            exon_len = exon.end - exon.start + 1

            counts = 0

            try:
                reads = bamfile.fetch(chrom, exon.start, exon.end)
            except ValueError:
                print "Error fetching region: %s:%d-%d" %(chrom,
                                                          exon.start,
                                                          exon.end)
                break
            
            # Count reads landing in exon
            for r in reads: counts += 1

            total_num_reads += counts

            # Skip exons that we've seen already or exons that are shorter
            # than the read length
            if (exon.start, exon.end) in regions_counted or \
               exon_len < read_len:
                continue
            
            exon_lens.append(exon_len)
            
            num_reads.append(counts)

            regions_counted[(exon.start, exon.end)] = True

        if len(regions_counted) == 0:
#            print "Gene %s exons are too small for %d-long reads" \
#                  %(gene_id, read_len)
            exons_too_small[gene_id] = total_num_reads
            continue

#        print "Used total of %d regions" %(len(regions_counted))
#        print "Total of %d regions are too small" %(num_too_small)

        rpkm = rpkm_per_region(exon_lens, num_reads, read_len,
                               num_total_reads)
        
#        print rpkm, exon_lens, num_reads, read_len      

        # Convert region lengths and number of reads to strings
        exon_lens_str = ",".join([str(e) for e in exon_lens])
        num_reads_str = ",".join([str(n) for n in num_reads])

        rpkm_entry = {'gene_id': gene_id,
                      'rpkm': "%.2f" %(rpkm),
                      'const_exon_lens': exon_lens_str,
                      'num_reads': num_reads_str}
        rpkms_dictlist.append(rpkm_entry)
        
#        print "RPKM: %.2f" %(rpkm)
            
        # Compute how many reads land in each constitutive exon
        num_genes += 1

    num_too_small = len(exons_too_small.keys())

    print "Computed RPKMs for %d genes." %(num_genes)
    print "  - Total of %d genes cannot be used because they lack const. regions." \
          %(num_no_const)
    print "  - Total of %d genes cannot be used since their exons are too small." \
          %(num_too_small)
    for gene, total_counts in exons_too_small.iteritems():
        print "      gene_id\ttotal_counts"
        print "    * %s\t%d" %(gene, total_counts)

    # Output RPKMs to file
    dictlist2file(rpkms_dictlist, output_filename,
                  rpkm_fieldnames)

    return rpkms_dictlist
    

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--compute-rpkm", dest="compute_rpkm", nargs=3, default=None,
                      help="Compute RPKMs.  Takes a GFF file with genes, an indexed/sorted BAM format "
                      "and an output directory.")
    parser.add_option("--read-len", dest="read_len", nargs=1, type="int", default=0,
                      help="Read length to use for RPKM computation.")
    (options, args) = parser.parse_args()

    if options.compute_rpkm != None:
        if options.read_len == 0:
            print "Error: Must give --read-len to compute RPKMs."
            return
        
        gff_filename = os.path.abspath(os.path.expanduser(options.compute_rpkm[0]))
        bam_filename = os.path.abspath(os.path.expanduser(options.compute_rpkm[1]))
        output_dir = os.path.abspath(os.path.expanduser(options.compute_rpkm[2]))
        
        compute_rpkm(gff_filename, bam_filename, options.read_len, output_dir)

if __name__ == "__main__":
    main()

