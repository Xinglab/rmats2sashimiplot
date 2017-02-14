# -*- mode: python; -*-
##
## Generate MISO summary files
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
import misopy.samples_utils as samples_utils

from misopy.parse_csv import *

import numpy as np
np.seterr(all='ignore')

miso_path = os.path.dirname(os.path.abspath(__file__))

def greeting(parser=None):
    print "MISO (Mixture of Isoforms model)"
    print "Summarize MISO output to get Psi values and confidence intervals."
    print "Use --help argument to view options.\n"
    if parser is not None:
        parser.print_help()
    
    
def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--summarize-samples", dest="summarize_samples",
                      nargs=2, default=None,
		      help="Compute summary statistics of the given set "
                      "of samples. Expects a directory with MISO output "
                      "and a directory to output summary file to.")
    parser.add_option("--summary-label", dest="summary_label",
                      nargs=1, default=None,
                      help="Label for MISO summary file. If not given, "
                      "uses basename of MISO output directory.")
    parser.add_option("--use-compressed", dest="use_compressed",
                      nargs=1, default=None,
                      help="Use compressed event IDs. Takes as input a "
                      "genes_to_filenames.shelve file produced by the "
                      "index_gff script.")
    (options, args) = parser.parse_args()

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
            misc_utils.make_dir(summary_output_dir)
	    
	summary_filename = os.path.join(summary_output_dir,
					'%s.miso_summary' %(samples_label))
	samples_utils.summarize_sampler_results(samples_dir,
                                            summary_filename,
                                            use_compressed=use_compressed)



if __name__ == "__main__":
    main()
