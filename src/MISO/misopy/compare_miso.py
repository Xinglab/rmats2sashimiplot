# -*- mode: python; -*-
##
## Comparing MISO outputs to get Delta Psi values and Bayes factors
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

def greeting(parser=None):
    print "MISO (Mixture of Isoforms model)"
    print "Compare MISO samples to get differential isoform statistics."
    print "Use --help argument to view options.\n"
    print "Example usage:\n"
    print "compare_miso --compare-samples sample1/ sample2/ results/"
    if parser is not None:
        parser.print_help()
    
    
def main():
    from optparse import OptionParser
    parser = OptionParser()
    
    ##
    ## Psi utilities
    ##
    parser.add_option("--compare-samples", dest="samples_to_compare",
                      nargs=3, default=None,
                      help="Compute comparison statistics between the two " \
                      "given samples. Expects three directories: the first is " \
                      "sample1's MISO output, the second is sample2's MISO " \
                      "output, and the third is the directory where " \
                      "results of the sample comparison will be outputted.")
    parser.add_option("--comparison-labels", dest="comparison_labels",
                      nargs=2, default=None,
                      help="Use these labels for the sample comparison "
                      "made by --compare-samples. "
                      "Takes two arguments: the label for sample 1 "
                      "and the label for sample 2, where sample 1 and "
                      "sample 2 correspond to the order of samples given "
                      "to --compare-samples.")
    parser.add_option("--use-compressed", dest="use_compressed",
                      nargs=1, default=None,
                      help="Use compressed event IDs. Takes as input a "
                      "genes_to_filenames.shelve file produced by the "
                      "index_gff script.")
    (options, args) = parser.parse_args()

    if options.samples_to_compare is None:
        greeting()

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


if __name__ == '__main__':
    main()

