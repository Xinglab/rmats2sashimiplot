##
## Draw gene structure from a GFF file
##
import os, sys, operator, subprocess
import math
import pysam
import glob
from pylab import *
from matplotlib.patches import PathPatch 
from matplotlib.path import Path

import misopy
import misopy.gff_utils as gff_utils
import misopy.sam_utils as sam_utils

from misopy.sashimi_plot.Sashimi import Sashimi
import misopy.sashimi_plot.plot_utils.plotting as plotting
import misopy.sashimi_plot.plot_utils.plot_settings as plot_settings
from misopy.sashimi_plot.plot_utils.plotting import show_spines
from misopy.parse_gene import parseGene

def plot_density_single(settings, sample_label,
                        tx_start, tx_end, gene_obj, mRNAs, strand,
                        graphcoords, graphToGene, bam_filename, axvar, chrom,
                        paired_end=False,
                        intron_scale=30,
                        exon_scale=4,
                        color='r',
                        ymax=None,
                        logged=False,
                        coverage=1,
                        number_junctions=True,
                        resolution=.5,
                        showXaxis=True,
                        showYaxis=True,
                        nyticks=3,
                        nxticks=4,
                        show_ylabel=True,
                        show_xlabel=True,
                        font_size=6,
                        junction_log_base=10,
                        plot_title=None,
                        plot_label=None):
    """
    Plot MISO events using BAM files and posterior distribution files.
    TODO: If comparison files are available, plot Bayes factors too.
    """
    bamfile = pysam.Samfile(bam_filename, 'rb')
    try:
        subset_reads = bamfile.fetch(reference=chrom, start=tx_start,end=tx_end)
    except ValueError as e:
        print "Error retrieving files from %s: %s" %(chrom, str(e))
        print "Are you sure %s appears in your BAM file?" %(chrom)
        print "Aborting plot..."
        return axvar
    wiggle, jxns = readsToWiggle_pysam(subset_reads, tx_start, tx_end)
    wiggle = 1e3 * wiggle / coverage
                
    # gene_reads = sam_utils.fetch_bam_reads_in_gene(bamfile, gene_obj.chrom,\
    #     tx_start, tx_end, gene_obj)
    # reads, num_raw_reads = sam_utils.sam_parse_reads(gene_reads,\
    #     paired_end=paired_end)
    # wiggle, jxns = readsToWiggle(reads, tx_start, tx_end)
    #wiggle = 1e3 * wiggle / coverage
 
    if logged:
        wiggle = log10(wiggle + 1)
    
    maxheight = max(wiggle)
    if ymax is None:
        ymax = 1.1 * maxheight
    else:
        ymax = ymax
    ymin = -.5 * ymax 

    # Reduce memory footprint by using incremented graphcoords.
    compressed_x = []
    compressed_wiggle = []
    prevx = graphcoords[0]
    tmpval = []
    for i in range(len(graphcoords)):
        tmpval.append(wiggle[i])
        if abs(graphcoords[i] - prevx) > resolution:
            compressed_wiggle.append(mean(tmpval))
            compressed_x.append(prevx)
            prevx = graphcoords[i]
            tmpval = []

    fill_between(compressed_x, compressed_wiggle,\
        y2=0, color=color, lw=0)
   
    sslists = []
    for mRNA in mRNAs:
        tmp = []
        for s, e in mRNA:
            tmp.extend([s, e])
        sslists.append(tmp)

    for jxn in jxns:
        leftss, rightss = map(int, jxn.split(":"))

        ss1, ss2 = [graphcoords[leftss - tx_start - 1],\
            graphcoords[rightss - tx_start]]

        mid = (ss1 + ss2) / 2
        h = -3 * ymin / 4
   
        numisoforms = 0
        for i in range(len(mRNAs)):
            if leftss in sslists[i] and \
                rightss in sslists[i]:
                numisoforms += 1
        if numisoforms > 0:
            if numisoforms % 2 == 0: # put on bottom 
                pts = [(ss1, 0), (ss1, -h), (ss2, -h), (ss2, 0)]
                midpt = cubic_bezier(pts, .5)
            else:                         # put on top 
                leftdens = wiggle[leftss - tx_start - 1]
                rightdens = wiggle[rightss - tx_start]

                pts = [(ss1, leftdens),
                       (ss1, leftdens + h),
                       (ss2, rightdens + h),
                       (ss2, rightdens)]
                midpt = cubic_bezier(pts, .5)

            if number_junctions:
                text(midpt[0], midpt[1], '%s'%(jxns[jxn]),
                     fontsize=6, ha='center', va='center', backgroundcolor='w')

            a = Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
            p = PathPatch(a, ec=color, lw=log(jxns[jxn] + 1) /\
                log(junction_log_base), fc='none')
            axvar.add_patch(p) 

    # Format plot
    # ylim(ymin, ymax)
    # axvar.spines['left'].set_bounds(0, ymax)
    axvar.spines['right'].set_color('none')
    axvar.spines['top'].set_color('none')

    if showXaxis:
        axvar.xaxis.set_ticks_position('bottom')
        xlabel('Genomic coordinate (%s), "%s" strand'%(gene_obj.chrom,
                                                       strand),
               fontsize=font_size)
        max_graphcoords = max(graphcoords) - 1
        coords_fontsize = font_size - (font_size * 0.2)
        xticks(linspace(0, max_graphcoords, nxticks),
               [graphToGene[int(x)] for x in \
                linspace(0, max_graphcoords, nxticks)],
               fontsize=coords_fontsize)
    else:
        axvar.spines['bottom'].set_color('none')
        xticks([])

    # if showYaxis:
    #     axvar.yaxis.set_ticks_position('left')
    #     yticks(linspace(0, ymax, nyticks), ['%d'%(x) for x in \
    #                                         linspace(0, ymax, nyticks)],
    #            fontsize=font_size)
    # else:
    #     axvar.spines['left'].set_color('none')
    #     yticks([])
    
    xlim(0, max(graphcoords))
    # Return modified axis
    return axvar


# Plot density for a series of bam files.
def plot_density(sashimi_obj, pickle_filename, event, plot_title=None):
#                 intron_scale=30, exon_scale=1, gene_posterior_ratio=5, posterior_bins=40,
#                 colors=None, ymax=None, logged=False, show_posteriors=True, coverages=None,
#                 number_junctions=True, resolution=.5, fig_width=8.5, fig_height=11,
#                 font_size=6, junction_log_base=10, reverse_minus=False,
#                 bar_posterior=False):
    # Get the settings we need
    settings = sashimi_obj.settings
    bam_files = settings["bam_files"]
    miso_files = settings["miso_files"]
    intron_scale = settings["intron_scale"]
    exon_scale = settings["exon_scale"]
    gene_posterior_ratio = settings["gene_posterior_ratio"]
    posterior_bins = settings["posterior_bins"]
    colors = settings["colors"]
    ymax = settings["ymax"]
    logged = settings["logged"]
    show_posteriors = settings["show_posteriors"]
    coverages = settings["coverages"]
    number_junctions = settings["number_junctions"]
    resolution = settings["resolution"]
    junction_log_base = settings["junction_log_base"]
    reverse_minus = settings["reverse_minus"]
    bar_posterior = settings["bar_posteriors"]
    font_size = settings["font_size"]
    nyticks = settings["nyticks"]
    nxticks = settings["nxticks"]
    show_ylabel = settings["show_ylabel"]
    show_xlabel = settings["show_xlabel"]
    if plot_title is None:
        plot_title = event
    print "Using intron scale ", intron_scale
    print "Using exon scale ", exon_scale

    # Always show y-axis for read densities for now
    showYaxis = True
    
    # Parse gene pickle file to get information about gene
    tx_start, tx_end, exon_starts, exon_ends, gene_obj, mRNAs, strand, chrom = \
        parseGene(pickle_filename, event)

    # Get the right scalings
    graphcoords, graphToGene = getScaling(tx_start, tx_end, strand,
                                          exon_starts, exon_ends, intron_scale,
                                          exon_scale, reverse_minus)

    nfiles = len(bam_files)
    if plot_title is not None:
        # Use custom title if given
        suptitle(plot_title, fontsize=10)
    else:
        suptitle(event, fontsize=10)
    plotted_axes = []
    
    for i in range(nfiles):
        if colors is not None:
            color = colors[i]
        else:
            color = None
        if coverages is not None:
            coverage = coverages[i]
        else:
            coverage = 1
        if i < nfiles - 1:
            showXaxis = False 
        else:
            showXaxis = True 

        bam_file = os.path.expanduser(bam_files[i])
        ax1 = subplot2grid((nfiles + 3, gene_posterior_ratio), (i, 0),
                           colspan=gene_posterior_ratio - 1)
        
        # Read sample label
        sample_label = settings["sample_labels"][i]

        print "Reading sample label: %s" %(sample_label)
        print "Processing BAM: %s" %(bam_file)
        
        plotted_ax = plot_density_single(settings, sample_label,
                                         tx_start, tx_end, gene_obj, mRNAs, strand,
                                         graphcoords, graphToGene, bam_file, ax1, chrom,
                                         paired_end=False, intron_scale=intron_scale,
                                         exon_scale=exon_scale, color=color,
                                         ymax=ymax, logged=logged, coverage=coverage,
                                         number_junctions=number_junctions, resolution=resolution,
                                         showXaxis=showXaxis, nyticks=nyticks, nxticks=nxticks,
                                         show_ylabel=show_ylabel, show_xlabel=show_xlabel,
                                         font_size=font_size,
                                         junction_log_base=junction_log_base)
        plotted_axes.append(plotted_ax)

        if show_posteriors:
            miso_file = os.path.expanduser(miso_files[i])
            try:
                ax2 = subplot2grid((nfiles + 3, gene_posterior_ratio),\
                    (i, gene_posterior_ratio - 1))

                if not os.path.isfile(miso_file):
                    print "Warning: MISO file %s not found" %(miso_file)

                print "Loading MISO file: %s" %(miso_file)
                plot_posterior_single(miso_file, ax2, posterior_bins,
                                      showXaxis=showXaxis, show_ylabel=False,
                                      font_size=font_size,
                                      bar_posterior=bar_posterior)
            except:
                box(on=False)
                xticks([])
                yticks([])
                print "Posterior plot failed."

    ##
    ## Figure out correct y-axis values
    ##
    ymax_vals = []
    if ymax != None:
        # Use user-given ymax values if provided
        max_used_yval = ymax
    else:
        # Compute best ymax value for all samples: take
        # maximum y across all.
        used_yvals = [curr_ax.get_ylim()[1] for curr_ax in plotted_axes]
        # Round up
        max_used_yval = math.ceil(max(used_yvals))

    # Reset axes based on this.
    # Set fake ymin bound to allow lower junctions to be visible
    fake_ymin = -0.6 * max_used_yval
    universal_yticks = linspace(0, max_used_yval,
                                nyticks + 1)
    # Round up yticks
    universal_ticks = map(math.ceil, universal_yticks)
    for sample_num, curr_ax in enumerate(plotted_axes):
        if showYaxis:
            curr_ax.set_ybound(lower=fake_ymin, upper=max_used_yval)
            curr_yticklabels = []
            for label in universal_yticks:
                if label <= 0:
                    # Exclude label for 0
                    curr_yticklabels.append("")
                else:
                    if label % 1 != 0:
                        curr_yticklabels.append("%.1f" %(label))
                    else:
                        curr_yticklabels.append("%d" %(label))
            curr_ax.set_yticklabels(curr_yticklabels,
                                    fontsize=font_size)
            curr_ax.spines["left"].set_bounds(0, max_used_yval)
            curr_ax.set_yticks(universal_yticks)
            curr_ax.yaxis.set_ticks_position('left')
            curr_ax.spines["right"].set_color('none')
            if show_ylabel:
                y_horz_alignment = 'left'
                if logged:
                    curr_ax.set_ylabel('RPKM $(\mathregular{\log}_{\mathregular{10}})$',
                                       fontsize=font_size,
                                       ha=y_horz_alignment)
                else:
                    curr_ax.set_ylabel('RPKM',
                                       fontsize=font_size,
                                       va="bottom",
                                       ha=y_horz_alignment)

        else:
            curr_ax.spines["left"].set_color('none')
            curr_ax.spines["right"].set_color('none')
            curr.ax.set_yticks([])
        ##
        ## Plot sample labels
        ##
        sample_color = colors[sample_num]
        # Make sample label y position be halfway between highest
        # and next to highest ytick
        if len(universal_yticks) >= 2:
            halfway_ypos = (universal_yticks[-1] - universal_yticks[-2]) / 2.
            label_ypos = universal_yticks[-2] + halfway_ypos
        else:
            label_ypos = universal_yticks[-1]
        curr_label = settings["sample_labels"][sample_num]
        curr_ax.text(max(graphcoords), label_ypos,
                     curr_label,
                     fontsize=font_size,
                     va='bottom',
                     ha='right',
                     color=sample_color)
                

    # Draw gene structure
    ax = subplot2grid((nfiles + 3, gene_posterior_ratio), (nfiles + 1, 0),
                      colspan=gene_posterior_ratio - 1, rowspan=2)
    plot_mRNAs(tx_start, mRNAs, strand, graphcoords, reverse_minus)
    subplots_adjust(hspace=.10, wspace=.7)


def getScaling(tx_start, tx_end, strand, exon_starts, exon_ends,
               intron_scale, exon_scale, reverse_minus):
    """
    Compute the scaling factor across various genic regions.
    """
   
    exoncoords = zeros((tx_end - tx_start + 1))
    for i in range(len(exon_starts)):
        exoncoords[exon_starts[i] - tx_start : exon_ends[i] - tx_start] = 1

    graphToGene = {}
    graphcoords = zeros((tx_end - tx_start + 1), dtype='f')
    x = 0
    if strand == '+' or not reverse_minus:
        for i in range(tx_end - tx_start + 1):
            graphcoords[i] = x
            graphToGene[int(x)] = i + tx_start
            if exoncoords[i] == 1:
                x += 1. / exon_scale
            else:
                x += 1. / intron_scale
    else:
        for i in range(tx_end - tx_start + 1):
            graphcoords[-(i + 1)] = x
            graphToGene[int(x)] = tx_end - i + 1
            if exoncoords[-(i + 1)] == 1:
                x += 1. / exon_scale
            else:
                x += 1. / intron_scale
    return graphcoords, graphToGene


def readsToWiggle_pysam(reads, tx_start, tx_end):
    """
    Convert reads to wiggles; uses pysam.
    """
    wiggle = zeros((tx_end - tx_start + 1), dtype='f')
    jxns = {}
    for read in reads:
        # Skip reads with no CIGAR string
        if read.cigar is None:
            print "Skipping read with no CIGAR string: %s" %(read.cigar)
            continue
        cigar_str = sam_utils.sam_cigar_to_str(read.cigar)

        if ("N" in cigar_str) and (cigar_str.count("N") > 1):
            print "Skipping read with multiple junctions crossed: %s" \
                  %(cigar_str)
            continue

        # Check if the read contains an insertion (I)
        # or deletion (D) -- if so, skip it
        for cigar_part in read.cigar:
            if cigar_part[0] == 1 or \
               cigar_part[1] == 2:
                print "Skipping read with CIGAR %s" \
                      %(cigar_str)
        aligned_positions = read.positions
        for i, pos in enumerate(aligned_positions):
            if pos < tx_start or pos > tx_end:
#                print "=>",pos
                continue
            wig_index = pos-tx_start
            wiggle[wig_index] += 1./read.qlen
            try:
                # if there is a junction coming up                
                if aligned_positions[i+1] > pos + 1: 
                    leftss = pos+1
                    rightss= aligned_positions[i+1]+1
                    if leftss > tx_start and leftss < tx_end \
                           and rightss > tx_start and rightss < tx_end:                      
                        jxn = ":".join(map(str, [leftss, rightss]))
                        try:
                            jxns[jxn] += 1 
                        except:
                            jxns[jxn] = 1
            except:
                pass
    return wiggle, jxns


# def readsToWiggle(reads, tx_start, tx_end):
#     """
#     Get wiggle and junction densities from reads.
#     """
#     read_positions, read_cigars = reads
#     wiggle = zeros((tx_end - tx_start + 1), dtype='f')
#     jxns = {}
#     for i in range(len(read_positions)):
#         pos, cigar = [read_positions[i], read_cigars[i]]
#         if "N" not in cigar:
#             rlen = int(cigar[:-1])
#             s = max([pos - tx_start, 0])
#             e = min([pos - tx_start + rlen, len(wiggle) - 1])
#             wiggle[s : e] += 1. / rlen
#         else:
#             left, right = cigar.split("N")
#             left, middle = map(int, left.split("M"))
#             right = int(right[:-1])
#             rlen = left + right
#             s1 = pos - tx_start
#             e1 = pos - tx_start + left
#             s2 = pos + left + middle - tx_start
#             e2 = pos + left + middle + right - tx_start

#             # Include read coverage from adjacent junctions.
#             if (e1 >= 0 and e1 < len(wiggle)) or (s1 >= 0 and s1 < len(wiggle)):
#                 wiggle[max([s1, 0]) : min([e1, len(wiggle)])] += 1. / rlen
#             if (e2 >= 0 and e2 < len(wiggle)) or (s2 >= 0 and s2 < len(wiggle)):
#                 wiggle[max([s2, 0]) : min([e2, len(wiggle)])] += 1. / rlen

#             # Plot a junction if both splice sites are within locus.
#             leftss = pos + left
#             rightss = pos + left + middle + 1
#             if leftss - tx_start >= 0 and leftss - tx_start < len(wiggle) \
#                 and rightss - tx_start >= 0 and rightss - tx_start < \
#                 len(wiggle): 

#                 jxn = ":".join(map(str, [leftss, rightss]))
#                 try:
#                     jxns[jxn] += 1 
#                 except:
#                     jxns[jxn] = 1 
#     return wiggle, jxns


def plot_mRNAs(tx_start, mRNAs, strand, graphcoords, reverse_minus):
    """
    Draw the gene structure.
    """
    yloc = 0 
    exonwidth = .3
    narrows = 50

    for mRNA in mRNAs:
        for s, e in mRNA:
            s = s - tx_start
            e = e - tx_start
            x = [graphcoords[s], graphcoords[e], graphcoords[e], graphcoords[s]]
            y = [yloc - exonwidth / 2, yloc - exonwidth / 2,\
                yloc + exonwidth / 2, yloc + exonwidth / 2]
            fill(x, y, 'k', lw=.5, zorder=20)

        # Draw intron.
        axhline(yloc, color='k', lw=.5)

        # Draw intron arrows.
        spread = .2 * max(graphcoords) / narrows
        for i in range(narrows):
            loc = float(i) * max(graphcoords) / narrows
            if strand == '+' or reverse_minus:
                x = [loc - spread, loc, loc - spread]
            else:
                x = [loc + spread, loc, loc + spread]
            y = [yloc - exonwidth / 5, yloc, yloc + exonwidth / 5]
            plot(x, y, lw=.5, color='k')

        yloc += 1 

    xlim(0, max(graphcoords)) 
    ylim(-.5, len(mRNAs) + .5)
    box(on=False)
    xticks([])
    yticks([]) 



def plot_posterior_single(miso_f, axvar, posterior_bins,
                          showXaxis=True,
                          showYaxis=True,
                          show_ylabel=True,
                          font_size=6,
                          bar_posterior=False):
    """
    Plot a posterior probability distribution for a MISO event.
    """
    posterior_bins = int(posterior_bins) 
    psis = [] 
    for line in open(miso_f):
        if not line.startswith("#") and not line.startswith("sampled"):
            psi, logodds = line.strip().split("\t")
            psis.append(float(psi.split(",")[0]))
  
    ci = .95 
    alpha = 1 - ci
    lidx = int(round((alpha / 2) * len(psis)) - 1)
    # the upper bound is the (1-alpha/2)*n nth smallest sample, where n is
    # the number of samples
    hidx = int(round((1 - alpha / 2) * len(psis)) - 1)
    psis.sort()
    clow, chigh = [psis[lidx], psis[hidx]]

    nyticks = 4

    if not bar_posterior:
        y, x, p = hist(psis, linspace(0, 1, posterior_bins),\
            normed=True, facecolor='k', edgecolor='w', lw=.2) 
        axvline(clow, ymin=.33, linestyle='--', dashes=(1, 1), color='#CCCCCC', lw=.5)
        axvline(chigh, ymin=.33, linestyle='--', dashes=(1, 1), color='#CCCCCC', lw=.5)
        axvline(mean(psis), ymin=.33, color='r')

        ymax = max(y) * 1.5
        ymin = -.5 * ymax
#             "$\Psi$ = %.2f\n$\Psi_{0.05}$ = %.2f\n$\Psi_{0.95}$ = %.2f" %\

        text(1, ymax,
             "$\Psi$ = %.2f\n[%.2f, %.2f]" % \
             (mean(psis), clow, chigh),
             fontsize=font_size,
             va='top',
             ha='left')

        ylim(ymin, ymax)
        axvar.spines['left'].set_bounds(0, ymax)
        axvar.spines['right'].set_color('none')
        axvar.spines['top'].set_color('none')
        axvar.spines['bottom'].set_position(('data', 0)) 
        axvar.xaxis.set_ticks_position('bottom')
        axvar.yaxis.set_ticks_position('left')
        if showYaxis:
            yticks(linspace(0, ymax, nyticks),\
                ["%d"%(y) for y in linspace(0, ymax, nyticks)],\
                fontsize=font_size)
        else:
            yticks([])
        if show_ylabel:
            ylabel("Frequency", fontsize=font_size, ha='right', va='center')
    else:
        ##
        ## Plot a horizontal bar version of the posterior distribution,
        ## showing only the mean and the confidence bounds.
        ##
        mean_psi_val = mean(psis)
        clow_err = mean_psi_val - clow
        chigh_err = chigh - mean_psi_val
        errorbar([mean_psi_val], [1],
                 xerr=[[clow_err], [chigh_err]],
                 fmt='o',
                 ms=4,
                 ecolor='k',
                 markerfacecolor="#ffffff",
                 markeredgecolor="k")
        text(1, 1,
             "$\Psi$ = %.2f\n[%.2f, %.2f]" % \
             (mean(psis), clow, chigh),
             fontsize=font_size,
             va='top',
             ha='left')
        yticks([])

    # Use same x-axis for all subplots
    # but only show x-axis labels for the bottom plot
    xlim([0, 1])
    psi_axis_fontsize = font_size - (font_size * 0.3)
    xticks([0, .2, .4, .6, .8, 1], fontsize=psi_axis_fontsize)
    
    if (not bar_posterior) and showYaxis:
        axes_to_show = ['bottom', 'left']
    else:
        axes_to_show = ['bottom']

    # Adjust x-axis to be lighter
    axis_size = 0.2
    tick_size = 1.2
    axis_color = "k"
    for shown_axis in axes_to_show:
        if shown_axis in axvar.spines:
            print "Setting color on %s axis" %(shown_axis)
            axvar.spines[shown_axis].set_linewidth(axis_size)
            axvar.xaxis.set_tick_params(size=tick_size,
                                        color=axis_color)
    if showXaxis:
        from matplotlib.ticker import FormatStrFormatter
        majorFormatter = FormatStrFormatter('%g')
        axvar.xaxis.set_major_formatter(majorFormatter)
        
        [label.set_visible(True) for label in axvar.get_xticklabels()]
        xlabel("MISO $\Psi$", fontsize=font_size)
        show_spines(axvar, axes_to_show)
    else:
        show_spines(axvar, axes_to_show)
        [label.set_visible(False) for label in axvar.get_xticklabels()]


def cubic_bezier(pts, t):
    """
    Get points in a cubic bezier.
    """
    p0, p1, p2, p3 = pts
    p0 = array(p0)
    p1 = array(p1)
    p2 = array(p2)
    p3 = array(p3)
    return p0 * (1 - t)**3 + 3 * t * p1 * (1 - t) ** 2 + \
        3 * t**2 * (1 - t) * p2 + t**3 * p3

    
def plot_density_from_file(settings_f, pickle_filename, event,
                           output_dir,
                           no_posteriors=False,
                           plot_title=None,
                           plot_label=None):
    """
    Read MISO estimates given an event name.
    """
    ##
    ## Read information about gene
    ##
    tx_start, tx_end, exon_starts, exon_ends, gene_obj, mRNAs, strand, chrom = \
        parseGene(pickle_filename, event)

    # Override settings flag on whether to show posterior plots
    # if --no-posteriors was given to plot.py
    sashimi_obj = Sashimi(event, output_dir,
                          event=event,
                          chrom=chrom,
                          settings_filename=settings_f,
                          no_posteriors=no_posteriors)
    
    print "Plotting read densities and MISO estimates along event..."
    print "  - Event: %s" %(event)
    
    settings = sashimi_obj.settings
    if no_posteriors:
        settings["show_posteriors"] = False
    bam_files = settings['bam_files']
    miso_files = settings['miso_files']

    # Setup the figure
    sashimi_obj.setup_figure()

    plot_density(sashimi_obj, pickle_filename, event,
                 plot_title=plot_title)

    # Save figure
    sashimi_obj.save_plot(plot_label=plot_label)
                 # intron_scale=settings["intron_scale"],
                 # exon_scale=settings["exon_scale"],
                 # gene_posterior_ratio=settings["gene_posterior_ratio"],
                 # posterior_bins=settings["posterior_bins"],
                 # show_posteriors=settings["show_posteriors"],
                 # logged=settings["logged"],
                 # colors=settings["colors"],
                 # ymax=settings["ymax"],
                 # coverages=settings["coverages"],
                 # number_junctions=settings["number_junctions"],
                 # resolution=settings["resolution"],
                 # fig_width=settings["fig_width"],
                 # fig_height=settings["fig_height"],
                 # font_size=settings["font_size"],
                 # junction_log_base=settings["junction_log_base"],
                 # reverse_minus=settings["reverse_minus"],
                 # bar_posterior=settings["bar_posteriors"])


