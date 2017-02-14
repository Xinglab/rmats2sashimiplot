#!/usr/bin/env python
import os
import sys
import unittest

import pysam
import sam_utils

class TestMISO(unittest.TestCase):
    """
    Test MISO functionality.
    """
    def setUp(self):
        # Find out the current directory
        self.miso_path = \
            os.path.dirname(os.path.abspath(os.path.expanduser(__file__)))
        self.tests_data_dir = \
            os.path.join(self.miso_path, "test-data")
        self.events_analysis_cmd = "miso"
        self.tests_output_dir = \
            os.path.join(self.miso_path, "test-output")
        self.test_sam_filename = \
            os.path.join(self.tests_data_dir,
                         "sam-data",
                         "c2c12.Atp2b1.sam")
        self.gff_events_dir = \
            os.path.join(self.miso_path, "gff-events")
        self.sam_to_bam_script = "sam_to_bam"
        self.index_gff_script = "index_gff"


    def test_a_sam_to_bam(self):
        """
        Test conversion of SAM to BAM.

        The 'a' ensures this runs first.
        """

        print "Testing conversion of SAM to BAM..."
        output_dir = \
            os.path.join(self.tests_output_dir, "sam-output")
        sam_to_bam_cmd = \
            "%s --convert %s %s" %(self.sam_to_bam_script,
                                   self.test_sam_filename,
                                   output_dir)
        print "Executing: %s" %(sam_to_bam_cmd)
        os.system(sam_to_bam_cmd)

        # Make sure conversion worked; sorted, indexed BAM file is outputted
        assert(os.path.exists(os.path.join(output_dir,
                                           "c2c12.Atp2b1.sorted.bam")))


    def test_a2_strandedness(self):
        """
        Test that strandedness is read correctly.
        """
        # Read 1 is forward, on plus strand
        # Has flag 129, i.e. '0b10000001'
        f_read = pysam.AlignedRead()
        f_read.qname = "f_read"
        f_read.flag = 129
        f_read.rname = 9
        f_read.pos = 4991443

        # Read 2 is reverse, on minus strand
        # Has flag 81, i.e. '0b1010001'
        r_read = pysam.AlignedRead()        
        r_read.qname = "r_read"
        r_read.flag = 81
        r_read.rname = 9
        r_read.pos = 4991578

        # Test that we can read the BAM strand flag correctly
        assert(sam_utils.flag_to_strand(f_read.flag) == "+"), \
            "Error in determining plus strand of read."
        assert(sam_utils.flag_to_strand(r_read.flag) == "-"), \
            "Error in determining minus strand of read."
        ##
        ## Test stranded-ness rules
        ##
        #   fr-unstranded,
        #   fr-firststrand,
        plus_target_strand = "+"
        minus_target_strand = "-"
        # fr-unstranded: both strand reads should match
        # either target strand
        print "Testing fr-unstranded..."
        for curr_read in [f_read, r_read]:
            for target in [plus_target_strand, minus_target_strand]:
                print "Checking read ", curr_read.qname, " against ", target
                assert(sam_utils.read_matches_strand(curr_read,
                                                     target,
                                                     "fr-unstranded") == True), \
                    "Error checking strand of fr-unstranded."
        # fr-firststrand: forward read must match target strand,
        # i.e. +read matches +target, and -read matches -target
        # test +read
        print "Testing fr-firststrand..."
        assert(sam_utils.read_matches_strand(f_read,
                                             plus_target_strand,
                                             "fr-firststrand") == True), \
            "+read must match +target under fr-firstrand."
        assert(sam_utils.read_matches_strand(f_read,
                                             minus_target_strand,
                                             "fr-firststrand") == False), \
            "+read must match +target under fr-firststrand."
        # test -read
        assert(sam_utils.read_matches_strand(r_read,
                                             plus_target_strand,
                                             "fr-firststrand") == False), \
            "-read must match -target under fr-firststrand."
        assert(sam_utils.read_matches_strand(r_read,
                                             minus_target_strand,
                                             "fr-firststrand") == True), \
            "-read must match -target under fr-firststrand."
        # Test fr-firststrand read pair
        pe = (300, 10)
        assert(sam_utils.read_matches_strand((f_read, r_read),
                                             plus_target_strand,
                                             "fr-firststrand",
                                             paired_end=pe) == True), \
            "(+, -) must match +target under fr-firststrand."
        # If target strand is -, second read must match minus strand
        assert(sam_utils.read_matches_strand((f_read, r_read),
                                             minus_target_strand,
                                             "fr-firststrand",
                                             paired_end=pe) == True), \
            "(+, -) must match -target under fr-firststrand."
        
    def test_z_gene_psi(self):
        """
        Test gene-level Psi inferences using SAM/BAM reads.

        The 'z' ensures this runs last.
        """
        print "Testing gene-level Psi..."
        sam_dir = os.path.join(self.tests_output_dir, "sam-output")
        bam_filename = os.path.join(sam_dir, "c2c12.Atp2b1.sorted.bam")

        read_len = 36
        insert_mean = 250
        insert_sd = 30

        # First index the GFF of interest
        gff_filename = os.path.join(self.gff_events_dir,
                                    "mm9",
                                    "genes",
                                    "Atp2b1.mm9.gff")
        gff_index_dir = os.path.join(self.gff_events_dir,
                                     "mm9",
                                     "genes",
                                     "Atp2b1",
                                     "indexed")
        print "Testing GFF indexing of: %s" %(gff_filename)
        index_cmd = "%s --index %s %s" %(self.index_gff_script,
                                         gff_filename,
                                         gff_index_dir)

        print "Executing: %s" %(index_cmd)
        os.system(index_cmd)

        output_dir = os.path.join(self.tests_output_dir,
                                  "gene-psi-output")
        miso_cmd = "%s --run %s %s --output-dir %s --read-len %d " \
                   %(self.events_analysis_cmd,
                     gff_index_dir,
                     bam_filename,
                     output_dir,
                     read_len)
        print "Executing: %s" %(miso_cmd)
        os.system(miso_cmd)

def main():
    unittest.main()
        
        
if __name__ == '__main__':
    main()
