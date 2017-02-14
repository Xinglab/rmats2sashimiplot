##
## Parse plotting configuration files for sashimi_plot
##
import sys
import os
import ast

import ConfigParser

import misopy
import misopy.miso_utils as miso_utils

def get_default_settings():
    """
    Return default settings for sashimi_plot.
    """
    settings = {"intron_scale": 30,
                "exon_scale": 1,
                "logged": False,
                "ymax": None,
                "show_posteriors": True,
                "number_junctions": True,
                "posterior_bins": 40,
                "gene_posterior_ratio": 5,
                "resolution": .5,
                "fig_width": 8.5,
                "fig_height": 11,
                "bar_posteriors": False,
                "junction_log_base": 10.,
                "reverse_minus": False,
                "bf_dist_bins": 20,
                "font_size": 6,
                "insert_len_bins": 25,
                "bf_thresholds": [0, 1, 2, 5, 10, 20],
                "nyticks": 3,
                "nxticks": 4,
                "show_ylabel": True,
                "show_xlabel": True,
                "sans_serif": False,
                "bar_color": "k"}
    return settings

def parse_plot_settings(settings_filename, event=None, chrom=None,
                        # Float parameters
                        FLOAT_PARAMS=["intron_scale",
                                      "exon_scale",
                                      "ymax",
                                      "resolution",
                                      "fig_width",
                                      "fig_height",
                                      "font_size",
                                      "junction_log_base"],
                        # Integer parameters
                        INT_PARAMS=["posterior_bins",
                                    "gene_posterior_ratio",
                                    "insert_len_bins",
                                    "nyticks",
                                    "nxticks"],
                        # Boolean parameters
                        BOOL_PARAMS=["logged",
                                     "show_posteriors",
                                     "number_junctions",
                                     "reverse_minus",
                                     "bar_posteriors",
                                     "show_ylabel",
                                     "show_xlabel",
                                     "sans_serif"],
                        # Parameters to be interpreted as Python lists or
                        # data structures
                        DATA_PARAMS=["miso_files",
                                     "bam_files",
                                     "bf_thresholds",
                                     "bar_color",
                                     "sample_labels"],
                        no_posteriors=False):
    """
    Populate a settings dictionary with the plotting parameters, parsed
    as the right datatype.
    """
    settings = get_default_settings()
    
    config = ConfigParser.ConfigParser()

    print "Reading settings from: %s" %(settings_filename)
    config.read(settings_filename)
    
    for section in config.sections():
        for option in config.options(section):
            print "Parsing %s:%s" %(section, option)
            if option in FLOAT_PARAMS:
                settings[option] = config.getfloat(section, option)
            elif option in INT_PARAMS:
                settings[option] = config.getint(section, option)
            elif option in BOOL_PARAMS:
                settings[option] = config.getboolean(section, option)
            elif option in DATA_PARAMS:
                settings[option] = ast.literal_eval(config.get(section,
                                                               option))
            else:
                settings[option] = config.get(section, option)

    # Ensure that bf_thresholds are integers
    settings["bf_thresholds"] = [int(t) for t in settings["bf_thresholds"]]
    
    if "colors" in settings:
        colors = ast.literal_eval(settings["colors"])
    else:
        colors = [None for x in settings["bam_files"]]
    settings["colors"] = colors
        
    if "bam_prefix" in settings:
        bam_files = [os.path.join(settings["bam_prefix"], x) \
                    for x in settings["bam_files"]]
    else:
        bam_files = settings["bam_files"]
    settings["bam_files"] = bam_files

    # Make the sample labels be the BAM files by default
    if "sample_labels" not in settings:
        settings["sample_labels"] = [os.path.basename(bfile) \
                                     for bfile in settings["bam_files"]]
    num_labels = len(settings["sample_labels"])
    num_bams = len(settings["bam_files"])
    num_colors = len(settings["colors"])
    if not (num_labels == num_bams == num_colors):
        print "Error: Must provide sample label and color for each entry in bam_files!"
        print "  - Provided %d labels, %d BAMs, %d colors" \
            %(num_labels, num_bams, num_colors)
        sys.exit(1)

    if no_posteriors:
        settings["show_posteriors"] = False

    if ("miso_prefix" in settings) and (event != None and chrom != None) \
        and settings["show_posteriors"]:
        miso_files = miso_utils.get_miso_output_files(event, chrom, settings)
    elif "miso_files" in settings:
        miso_files = settings["miso_files"]
    else:
        miso_files = []
    settings["miso_files"] = miso_files
    
    if "coverages" in settings:
        coverages = ast.literal_eval(settings["coverages"])
        coverages = map(float, coverages)
        # Normalize coverages per M
        coverages = [x / 1e6  for x in coverages]
    else:
        coverages = [1 for x in settings["bam_files"]]
    settings["coverages"] = coverages

    if len(settings["coverages"]) != len(settings["sample_labels"]):
        print "Error: Must provide a coverage value for each sample or leave coverages unset."
        sys.exit(1)
    
    return settings
