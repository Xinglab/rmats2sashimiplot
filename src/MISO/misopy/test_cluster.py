#!/usr/bin/env python
##
## Test running of MISO on cluster
##
import unittest
import os

class TestCluster(unittest.TestCase):
    """
    Test cluster functionality.
    """
    def setUp(self):
        # Find out the current directory
        self.miso_path = os.path.dirname(os.path.abspath(os.path.expanduser(__file__)))
        self.tests_data_dir = os.path.join(self.miso_path, "test-data")
        self.events_analysis_cmd = "python %s " %(os.path.join(self.miso_path,
                                                               "run_events_analysis.py"))
        self.tests_output_dir = os.path.join(self.miso_path, "test-output")
        self.test_sam_filename = os.path.join(self.tests_data_dir,
                                              "sam-data",
                                              "c2c12.Atp2b1.sam")
        self.gff_events_dir = os.path.join(self.miso_path, "gff-events")
        self.sam_to_bam_script = os.path.join(self.miso_path, "sam_to_bam.py")
        self.index_gff_script = os.path.join(self.miso_path, "index_gff.py")

    def test_cluster_single_end_run(self):
        """
        Test MISO on cluster.
        """
        print "Testing single-end SE event interface..."

        ##
        ## Try running MISO on cluster using default settings.
        ##
        sample_name = "se-sample"
        counts_filename = os.path.join(self.tests_data_dir,
                                       "se-counts",
                                       "se_test.counts")
        output_dir = os.path.join(self.tests_output_dir, "SE-output")

        read_len = 35
        overhang_len = 4

        event_type = "SE"
        
        miso_cmd = "%s --compute-events-psi %s %s --output-dir %s --read-len %d --overhang-len %d " \
                   " --event-type %s --use-cluster " %(self.events_analysis_cmd,
                                                       sample_name,
                                                       counts_filename,
                                                       output_dir,
                                                       read_len,
                                                       overhang_len,
                                                       event_type)
        print "Executing: %s" %(miso_cmd)
        os.system(miso_cmd)

    def test_cluster_gene_psi(self):
        """
        Test gene-level Psi inferences using SAM/BAM reads.
        """
        print "Testing gene-level Psi..."
        sam_dir = os.path.join(self.tests_output_dir, "sam-output")
        bam_filename = os.path.join(sam_dir, "c2c12.Atp2b1.sorted.bam")

        read_len = 36
        insert_mean = 250
        insert_sd = 30

        # First index the GFF of interest
        gff_filename = os.path.join(self.gff_events_dir, "mm9", "genes", "Atp2b1.mm9.gff")
        gff_index_dir = os.path.join(self.gff_events_dir, "mm9", "indexed")
        index_cmd = "python %s --index %s %s" %(self.index_gff_script,
                                                gff_filename,
                                                gff_index_dir)

        print "Executing: %s" %(index_cmd)
        os.system(index_cmd)

        output_dir = os.path.join(self.tests_output_dir, "gene-psi-output")
        
        miso_cmd = "%s --compute-genes-psi %s %s --output-dir %s --read-len %d " \
                   " --paired-end %d %d --use-cluster" \
                   %(self.events_analysis_cmd,
                     gff_index_dir,
                     bam_filename,
                     output_dir,
                     read_len,
                     insert_mean,
                     insert_sd)
        print "Executing: %s" %(miso_cmd)
        os.system(miso_cmd)

        
if __name__ == '__main__':
    unittest.main()
