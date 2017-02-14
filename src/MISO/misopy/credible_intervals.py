from scipy import *
from numpy import *

def format_credible_intervals(event_name, samples,
                              confidence_level=0.95):
    """
    Returns a list of print-able credible intervals for an NxM samples
    matrix. Handles both the two isoform and multi-isoform cases.
    """
    num_samples, num_isoforms = shape(samples)
    if num_isoforms > 2:
        cred_interval = compute_multi_iso_credible_intervals(samples,
                                                             confidence_level=confidence_level)
        cred_interval_lowbounds = ",".join([str("%.2f" %(ci[0])) for ci in cred_interval])
        cred_interval_highbounds = ",".join([str("%.2f" %(ci[1])) for ci in cred_interval])
        posterior_mean = ",".join("%.2f" %(val) for val in mean(samples, 0))
        output_fields = [event_name,
                         "%s" %(posterior_mean),
                         "%s" %(cred_interval_lowbounds),
                         "%s" %(cred_interval_highbounds)]
    else:
        cred_interval = compute_credible_intervals(samples)
        posterior_mean = mean(samples, 0)[0]
        output_fields = [event_name,
                         "%.2f" %(posterior_mean),
                         "%.2f" %(cred_interval[0]),
                         "%.2f" %(cred_interval[1])]
    return output_fields


def compute_credible_intervals(samples, confidence_level=.95):
    """
    Compute Bayesian confidence intevals (credible intervals) for the set of samples given
    based on the method of Chen and Shao (1998).

    Assumes that samples is an Nx2 vector of posterior samples.
    """
    if samples.ndim == 2:
	samples = samples[:, 0]
    num_samples = len(samples)
    # confidence percentage is 100(1-alpha)%    
    alpha = 1 - confidence_level
    # compute the lower bound of the interval
    # the lower bound is the (alpha/2)*n-th smallest sample, where n is the
    # number of samples
    lower_bound_indx = round((alpha/2)*num_samples) - 1
    # the upper bound is the (1-alpha/2)*n nth smallest sample, where n is
    # the number of samples
    upper_bound_indx = round((1-alpha/2)*num_samples) - 1
    assert(lower_bound_indx > 0)
    assert(upper_bound_indx > 0)
    # sort samples along first axis
    samples.sort()
    cred_interval = [samples[lower_bound_indx], samples[upper_bound_indx]]
    return cred_interval


def compute_multi_iso_credible_intervals(multi_iso_samples,
                                         confidence_level=0.95):
    """
    Compute multiple isoforms credible intervals for a set of NxM matrix.
    """
    credible_intervals = []
    num_samples, num_isoforms = shape(multi_iso_samples)

    for iso_num in range(num_isoforms):
        ci = compute_credible_intervals(multi_iso_samples[:, iso_num],
                                        confidence_level=confidence_level)
        credible_intervals.append(ci)
        
    return credible_intervals
        
