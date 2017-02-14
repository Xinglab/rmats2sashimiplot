##
## Visualize samples produced by MISO.
##
## TODO: In future interface with spliceplot to produce densities along a gene model
##

from scipy import *
from numpy import *

import matplotlib
#from plotting import colors, show_spines, axes_square

import matplotlib.pyplot as plt
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
import time

import misopy
from misopy.parse_csv import csv2array
from misopy.sashimi_plot.plot_utils.plotting import *
import misopy.Gene as Gene
import misopy.hypothesis_test as ht

class SamplesPlotter:
    """
    Visualize a set of samples from a run of MISO.
    """
    def __init__(self, samples, params, log_scores=None,
                 percent_acceptance=None,
		 true_psi=None):
	"""
	Given a sampler instance, store its properties.
	"""
	# sampler parameters
	self.samples = samples	
	self.params = params
	self.log_scores = log_scores
	self.percent_acceptance = percent_acceptance
	self.true_psi = true_psi
	
	assert(len(samples) > 1)
        

    def plot(self, fig=None, output_dir=None, num_rows=1, num_cols=1, subplot_start=1,
             title=None, plot_intervals=None, value_to_label=None, label=None, bins=10,
             bbox_coords=None, vanilla=False,
             plot_mean=False, fig_dims=(5, 5)):
	"""
	Plot a set of samples.

	 - credible_intervals: if set to true, plot Bayesian confidence intervals 
	"""
	plot_handle = None

        num_samples, num_isoforms = shape(self.samples)

        if num_isoforms == 2:
            plot_handle = self.plot_two_iso_samples(fig=fig, plots_dir=output_dir, num_cols=num_cols,
                                                    num_rows=num_rows, subplot_start=subplot_start,
                                                    plot_intervals=plot_intervals,
                                                    value_to_label=value_to_label,
                                                    label=label, bbox_coords=bbox_coords,
                                                    title=title, vanilla=vanilla,
                                                    plot_mean=plot_mean, fig_dims=fig_dims)
	elif num_isoforms > 2:
	    num_isoforms = self.samples.shape[1] 
	    num_rows = 1
	    num_cols = num_isoforms
	    for c in range(num_cols):
		plot_handle = self.plot_two_iso_samples(fig, isoform_index=c,
                                                        subplot_start=c + 1, num_cols=num_cols,
                                                        plot_intervals=plot_intervals,
							title=title, bins=bins, vanilla=vanilla,
                                                        plot_mean=plot_mean, fig_dims=fig_dims)
		plt.ylabel('Frequency (Isoform %d)' %(c + 1))
	    plt.subplots_adjust(wspace=0.5)
        else:
            raise Exception, "Invalid number of isoforms %d" %(num_isoforms)
	return plot_handle
	
    def plot_two_iso_samples(self, fig=None, isoform_index=0, num_rows=1, num_cols=1, subplot_start=1,
			     plots_dir=None, map_estimate=None, simulation_num=1,
			     plot_intervals=False, value_to_label=None, label=None, plot_filename=None,
                             bins=None, bbox_coords=None, with_legend=True, title=None, vanilla=False,
                             plot_mean=False, normed=False, fig_dims=(5, 5)):
	"""
	Plot a set of samples for Psi of a two isoform gene.
	"""
	if not fig:
	    sampled_psi_fig = plt.figure(figsize=fig_dims, dpi=300)
	else:
	    sampled_psi_fig = fig
	ax = sampled_psi_fig.add_subplot(num_rows, num_cols, subplot_start)
	num_iters = int(self.params['iters'])
	burn_in = int(self.params['burn_in'])
	lag = int(self.params['lag'])
	percent_acceptance = float(self.params['percent_accept'])
	proposal_type = self.params['proposal_type']
	plt.rcParams['font.size'] = 10
	show_spines(ax, ['left', 'bottom'])
	bins = bins
	assert((value_to_label == None and label == None) or \
	       (value_to_label != None and label != None))
	# retrieve samples
	samples_to_plot = self.samples[:, isoform_index]
	# picasso blue #0276FD
        
	if not vanilla:
	    if bins != None:
		plt.hist(samples_to_plot, align='mid', lw=0.5, facecolor='#0276FD',
                         edgecolor='#ffffff')
	    else:
		plt.hist(samples_to_plot, align='mid', lw=0.5, facecolor='#0276FD',
                         edgecolor='#ffffff')
	else:
	    plt.hist(samples_to_plot, align='mid', facecolor='#0276FD', edgecolor='#0276FD')
	plt.xlabel(r'${\hat{\Psi}}_{\mathregular{MISO}}$')
	plt.ylabel('Frequency')
	plt.xlim([0, 1])

        # Normalize samples
        if normed:
            yticks = list(plt.gca().get_yticks())
            print "yticks: ", yticks
            ytick_labels = ["%.2f" %(float(ytick) / float(normed)) for ytick in yticks]
            ax.set_yticklabels(ytick_labels)
#            samples_to_plot = samples_to_plot / float(len(samples_to_plot))
        

            # curr_tick_labels = [label.get_label() for label in ax.get_yticklabels()]
            # print "Current tick labels: ", curr_tick_labels
            # new_tick_labels = []
            # for label in curr_tick_labels:
            #     if len(label) > 0:
            #         new_label = "%.1f" %(float(label) / normed)
            #     else:
            #         new_label = ""
            #     new_tick_labels.append(new_label)
            # #ax.set_yticklabels(new_tick_labels)
        
	curr_axes = plt.gca()
	# Plot MAP estimate for same data
	if map_estimate:
	    l = plt.axvline(x=map_estimate, color='b', linewidth=1.2, ls='-', label=r'${\hat{\Psi}}_{MAP}\ =\ %.2f$' %(map_estimate))
	# Plot true Psi
	if self.true_psi:
	    plot_id = "%dsimul_%diters_%dburnin_%dlag_%s_truepsi_%.2f.pdf" \
                      %(simulation_num, num_iters, burn_in, lag, proposal_type, self.true_psi)
	    l = plt.axvline(x=self.true_psi, color='r', linewidth=1.2, ls='-', label=r'True $\Psi$')
	else:
	    # Unknown true Psi
	    plot_id = "%dsimul_%diters_%dburnin_%dlag_%s_%s_truepsi.pdf" \
                      %(simulation_num, num_iters, burn_in, lag, proposal_type, 'unknown')
	if value_to_label:
	    l = plt.axvline(x=value_to_label, color='r', linewidth=1.2, ls='-', label=label)
	# plot credible intervals if given
	if plot_intervals:
#	    print "Plotting %.2f confidence intervals" %(plot_intervals * 100)
	    interval_c1, interval_c2 = ht.compute_credible_intervals(samples_to_plot, plot_intervals)
	    plt.axvline(x=interval_c1, color='#999999', linewidth=0.7, ls='--',
			label=r'%d' %(plot_intervals*100) + '% CI')
	    plt.axvline(x=interval_c2, color='#999999', linewidth=0.7, ls='--')
	if plot_mean:
	    sample_mean = mean(samples_to_plot)
	    plt.axvline(x=sample_mean, color='r', linewidth=0.8, label='Mean')
	if with_legend and (plot_intervals or self.true_psi):
	    if not bbox_coords:
		lg = plt.legend(handletextpad=0.172, borderpad=0.01, labelspacing=.008,
                                handlelength=1.4, loc='best', numpoints=1)
	    else:
		lg = plt.legend(handletextpad=0.172, borderpad=0.01, labelspacing=.008,
                                handlelength=1.4, loc='best', numpoints=1,
				bbox_to_anchor=bbox_coords)
	    lg.get_frame().set_linewidth(0)
	    for t in lg.get_texts():
		t.set_fontsize(8)
	if title:
	    plt.title(title)
	if plots_dir:
	    if not plot_filename:
		plt.savefig(plots_dir + "sampled_psi_hist_%s" %(plot_id))
	    else:
		plt.savefig(plots_dir + plot_filename + '.pdf')
	return curr_axes
	# Plot joint scores as function of number of samples
	#log_joint_fig = plt.figure(figsize=(7,4.5), dpi=300)
	#skip = 15
	#print "Skip of %d when plotting log joint scores" %(skip)
	#plt.plot(arange(0, len(total_log_scores), skip),
	#	 total_log_scores[arange(0, len(total_log_scores), skip)])
	#print "Total log scores plotted: ", len(total_log_scores)
	#plt.xlabel('Number of iterations (lag not shown)')
	#plt.ylabel('Log joint score')
	#plt.savefig(plots_dir + "log_joint_scores_skip%d_%s" %(skip, plot_id))
	    
# def load_samples(samples_filename):
#     """
#     Load a set of samples from a file and build an associated gene.
#     """
#     samples_data, h = csv2array(samples_filename, skiprows=1,
#                                 raw_header=True)
#     samples = []
#     log_scores = []
#     for line in samples_data:
# 	psi_vals = [float(v) for v in line['sampled_psi'].split(',')]
# 	samples.append(psi_vals)
# 	log_scores.append(float(line['log_score']))
#     params, gene = parse_sampler_params(h[0])
#     return (array(samples), array(log_scores), params, gene)

