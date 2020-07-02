# -*- mode: python; -*-
##
## sashimi_plot
##
## Utility for visualizing RNA-Seq densities along gene models and
## for plotting MISO output
##
import os
import sys
import glob
import matplotlib

# Add misopy path
miso_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, miso_path)

# Use PDF backend
matplotlib.use("pdf")

from scipy import *
from numpy import *

import pysam
import shelve

import misopy
import misopy.gff_utils as gff_utils
import misopy.pe_utils as pe_utils
from misopy.parse_csv import csv2dictlist_raw

from misopy.samples_utils import load_samples
from misopy.sashimi_plot.Sashimi import Sashimi
from misopy.sashimi_plot.plot_utils.samples_plotter import SamplesPlotter
from misopy.sashimi_plot.plot_utils.plotting import *
from misopy.sashimi_plot.plot_utils.plot_gene import plot_density_from_file
import matplotlib.pyplot as plt
from matplotlib import rc


def plot_bf_dist(bf_filename, settings_filename, output_dir,
                 max_bf=1e12):
    """
    Plot a Bayes factor distribution from a .miso_bf file.
    """
    if not bf_filename.endswith(".miso_bf"):
        print "WARNING: %s does not end in .miso_bf, are you sure it is the " \
              "output of a MISO samples comparison?" %(bf_filename)

    # Load BF data
    data, h = csv2dictlist_raw(bf_filename)

    plot_name = os.path.basename(bf_filename)
    sashimi_obj = Sashimi(plot_name, output_dir,
                          settings_filename=settings_filename)
    settings = sashimi_obj.settings

    # Setup the figure
    sashimi_obj.setup_figure()

    # Matrix of bayes factors and delta psi pairs
    bfs_and_deltas = []
    for event in data:
        bf = event['bayes_factor']
        delta_psi = event['diff']

        if type(bf) == str and "," in bf:
            print "WARNING: %s is a multi-isoform event, skipping..." \
                %(event)
            continue
        else:
            # Impose upper limit on Bayes factor
            bf = min(1e12, float(bf))
            delta_psi = float(delta_psi)

        bfs_and_deltas.append([bf, delta_psi])

    bfs_and_deltas = array(bfs_and_deltas)
    num_events = len(bfs_and_deltas)

    print "Loaded %d event comparisons." %(num_events)

    output_filename = sashimi_obj.output_filename

    print "Plotting Bayes factors distribution"
    print "  - Output filename: %s" %(output_filename)
    bf_thresholds = settings["bf_thresholds"]
    bar_color = settings["bar_color"]

    min_bf_thresh = min(bf_thresholds)
    num_events_used = sum(bfs_and_deltas[:, 0] >= min_bf_thresh)
    for thresh in bf_thresholds:
        if type(thresh) != int:
            print "Error: BF thresholds must be integers."
            sys.exit(1)
    print "Using BF thresholds: "
    print bf_thresholds
    print "Using bar color: %s" %(bar_color)
    plot_cumulative_bars(bfs_and_deltas[:, 0],
                         bf_thresholds,
                         bar_color=bar_color,
                         logged=True)
    plt.xticks(bf_thresholds)
    c = 1
    plt.xlim([bf_thresholds[0] - c, bf_thresholds[-1] + c])
    plt.title("Bayes factor distributions\n(using %d/%d events)" \
              %(num_events_used, num_events))
    plt.xlabel("Bayes factor thresh.")
    plt.ylabel("No. events")
    sashimi_obj.save_plot()



def plot_event(event_name, pickle_dir, settings_filename,
               output_dir,
               group_info=None,
               no_posteriors=False,
               plot_title=None,
               plot_label=None):

    """
    Visualize read densities across the exons and junctions
    of a given MISO alternative RNA processing event.

    Also plots MISO estimates and Psi values.
    """
    if not os.path.isfile(settings_filename):
        print "Error: settings filename %s not found." %(settings_filename)
        sys.exit(1)

    if not os.path.isdir(pickle_dir):
        print "Error: event pickle directory %s not found." %(pickle_dir)
        sys.exit(1)

    # Retrieve the full pickle filename
    genes_filename = os.path.join(pickle_dir,
                                  "genes_to_filenames.shelve")

    # Check that file basename exists
    if len(glob.glob("%s*" %(genes_filename))) == 0:
        raise Exception, "Cannot find file %s. Are you sure the events " \
                         "were indexed with the latest version of index_gff.py?" \
                         %(genes_filename)

    event_to_filenames = shelve.open(genes_filename)
    if event_name not in event_to_filenames:
        raise Exception, "Event %s not found in pickled directory %s. " \
              "Are you sure this is the right directory for the event?" \
              %(event_name, pickle_dir)

    pickle_filename = event_to_filenames[event_name]

    if no_posteriors:
        print "Asked to not plot MISO posteriors."

    plot_density_from_file(settings_filename, pickle_filename, event_name,
                           output_dir,
                           group_info=group_info,
                           no_posteriors=no_posteriors,
                           plot_title=plot_title,
                           plot_label=plot_label)


def plot_insert_len(insert_len_filename,
                    settings_filename,
                    output_dir):
    """
    Plot insert length distribution.
    """
    if not os.path.isfile(settings_filename):
        print "Error: settings filename %s not found." %(settings_filename)
        sys.exit(1)
    plot_name = os.path.basename(insert_len_filename)
    sashimi_obj = Sashimi(plot_name, output_dir,
                          settings_filename=settings_filename)
    settings = sashimi_obj.settings
    num_bins = settings["insert_len_bins"]
    output_filename = sashimi_obj.output_filename
    sashimi_obj.setup_figure()
    s = plt.subplot(1, 1, 1)

    print "Plotting insert length distribution..."
    print "  - Distribution file: %s" %(insert_len_filename)
    print "  - Output plot: %s" %(output_filename)

    insert_dist, params = pe_utils.load_insert_len(insert_len_filename)

    mean, sdev, dispersion, num_pairs \
          = pe_utils.compute_insert_len_stats(insert_dist)
    print "min insert: %.1f" %(min(insert_dist))
    print "max insert: %.1f" %(max(insert_dist))
    plt.title("%s (%d read-pairs)" \
              %(plot_name,
                num_pairs),
              fontsize=10)
    plt.hist(insert_dist, bins=num_bins, color='k',
             edgecolor="#ffffff", align='mid')
    axes_square(s)
    ymin, ymax = s.get_ylim()
    plt.text(0.05, 0.95, "$\mu$: %.1f\n$\sigma$: %.1f\n$d$: %.1f" \
             %(round(mean, 2),
               round(sdev, 2),
               round(dispersion, 2)),
             horizontalalignment='left',
             verticalalignment='top',
             bbox=dict(edgecolor='k', facecolor="#ffffff",
                       alpha=0.5),
             fontsize=10,
             transform=s.transAxes)
    plt.xlabel("Insert length (nt)")
    plt.ylabel("No. read pairs")
    sashimi_obj.save_plot()

def greeting():
    print "Sashimi plot: Visualize spliced RNA-Seq reads along gene models. " \
          "Part of the MISO (Mixture of Isoforms model) framework."
    print "See --help for usage.\n"
    print "Manual available at: http://genes.mit.edu/burgelab/miso/docs/sashimi.html\n"


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--plot-insert-len", dest="plot_insert_len", nargs=2, default=None,
                      help="Plot the insert length distribution from a given insert length (*.insert_len) "
                      "filename. Second argument is a settings file name.")
    parser.add_option("--plot-bf-dist", dest="plot_bf_dist", nargs=2, default=None,
                      help="Plot Bayes factor distributon. Takes the arguments: "
                      "(1) Bayes factor filename (*.miso_bf) filename, "
                      "(2) a settings filename.")
    parser.add_option("--plot-event", dest="plot_event", nargs=3, default=None,
                      help="Plot read densities and MISO inferences for a given alternative event. "
                      "Takes the arguments: (1) event name (i.e. the ID= of the event based on MISO gff3 "
                      "annotation file, (2) directory where indexed GFF annotation is (output of "
                      "index_gff.py), (3) path to plotting settings file.")
    parser.add_option("--no-posteriors", dest="no_posteriors", default=False, action="store_true",
                      help="If given this argument, MISO posterior estimates are not plotted.")
    parser.add_option("--plot-title", dest="plot_title", default=None, nargs=1,
                      help="Title of plot: a string that will be displayed at top of plot. Example: " \
                      "--plot-title \"My favorite gene\".")
    parser.add_option("--plot-label", dest="plot_label", default=None, nargs=1,
                      help="Plot label. If given, plot will be saved in the output directory as " \
                      "the plot label ending in the relevant extension, e.g. <plot_label>.pdf. " \
                      "Example: --plot-label my_gene")
    parser.add_option("--output-dir", dest="output_dir", nargs=1, default=None,
                      help="Output directory.")
    parser.add_option("--group-info", dest="group_info", nargs=1, default= None,
                      help="If there is the need to divide bam files into groups, then provided this parameter with the"
                           " the group files' name. Exemple:"
                           " \'--group-info gf.gf\'")  # TODO: modify it when the format is determined
    (options, args) = parser.parse_args()

    if options.plot_event is None:
        greeting()
        sys.exit(1)

    if options.output_dir == None:
        print "Error: need --output-dir"
        sys.exit(1)

    output_dir = os.path.abspath(os.path.expanduser(options.output_dir))

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    no_posteriors = options.no_posteriors

    plot_title = options.plot_title
    plot_label = options.plot_label

    if options.plot_insert_len != None:
        insert_len_filename = os.path.abspath(os.path.expanduser(options.plot_insert_len[0]))
        settings_filename = os.path.abspath(os.path.expanduser(options.plot_insert_len[1]))
        plot_insert_len(insert_len_filename, settings_filename, output_dir)

    if options.plot_bf_dist != None:
        bf_filename = os.path.abspath(os.path.expanduser(options.plot_bf_dist[0]))
        settings_filename = os.path.abspath(os.path.expanduser(options.plot_bf_dist[1]))
        plot_bf_dist(bf_filename, settings_filename, output_dir)

    if options.plot_event != None:
        event_name = options.plot_event[0]
        pickle_dir = os.path.abspath(os.path.expanduser(options.plot_event[1]))
        settings_filename = os.path.abspath(os.path.expanduser(options.plot_event[2]))
        group_info = options.group_info
        plot_event(event_name, pickle_dir, settings_filename, output_dir,
                   group_info=group_info,
                   no_posteriors=no_posteriors,
                   plot_title=plot_title,
                   plot_label=plot_label)


if __name__ == '__main__':
    main()
