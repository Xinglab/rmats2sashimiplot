# -*- mode: python; -*-
##
## Check if all the necessary modules to run MISO are available
##
import time
import os

import misopy
import misopy.misc_utils as utils

def check_module_availability(required_modules):
    unavailable_mods = 0
    print "Checking availability of Python modules for MISO"
    print "Looking for required Python modules.."
    for module_name in required_modules:
	print "Checking for availability of: %s" %(module_name)
	try:
	    __import__(module_name)
            # Manually check for correct matplotlib version
            # required for sashimi_plot
            if module_name == "matplotlib":
                import matplotlib.pyplot as plt
                if not hasattr(plt, "subplot2grid"):
                    print "WARNING: subplot2grid function is not available in matplotlib. " \
                          "to use sashimi_plot, you must upgrade your matplotlib " \
                          "to version 1.1.0 or later. This function is *not* required " \
                          "for MISO use."
	except ImportError:
	    print "  - Module %s not available!" %(module_name)
            if module_name == "matplotlib":
                print "matplotlib is required for sashimi_plot"
	    unavailable_mods += 1
    if unavailable_mods != 0:
        print "Total of %d modules were not available. " \
              "Please install these and try again." %(unavailable_mods)
    else:
        print "All modules are available!"
    print "Looking for required executables.."
    required_programs = ["samtools", "bedtools"]
    for prog in required_programs:
        p = utils.which(prog)
        print "Checking if %s is available" %(prog)
        if p is None:
            print " - Cannot find %s!" %(prog)
            if prog == "bedtools":
                print "bedtools is only required for prefiltering " \
                      "and computation of insert lengths."
                if utils.which("tagBam"):
                    print "Your bedtools installation might be available " \
                          "but outdated. Please upgrade bedtools and " \
                          "ensure that \'bedtools\' is available on path."
        else:
            print "  - %s is available" %(prog)
    return unavailable_mods


def main():
    required_modules = ['numpy', 'scipy', 'json', 'matplotlib',
                        'pysam']
    check_module_availability(required_modules)
    

if __name__ == '__main__':
    main()
    
