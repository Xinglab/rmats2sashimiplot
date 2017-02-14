# -*- mode: python; -*-
##
## Filter two-isoform events
##
import os
import sys
import time
import re
import string
import misopy
from misopy.parse_csv import *
from collections import defaultdict

def fname_callback(option, op_str, value, parser):
    """
    Handles the parsing the variable number of filenames to the --filter
    and --control flags, returns a list of the filenames
    """
    flist = []
    for arg in parser.rargs:
        if arg[0] is "-":
            break
        else:
            flist.append(arg)
    setattr(parser.values, option.dest, flist)

def get_counts(counts_str):
    """
    Return number of inclusion, exclusion and reads supporting
    both isoforms.
    """
    num_inc = 0
    num_exc = 0
    num_both = 0

    fields = re.findall("(\(.{3}\):\d+)", counts_str)

    isoforms = re.findall("\([01,]+\)", counts_str)[0]
    isoforms = isoforms.translate(None, string.punctuation)

    if len(isoforms) > 2:
        return None

    if len(fields) == 0:
        return None
    
    for field in fields:
        iso_type, count = field.split(":")
        count = int(count)

        if iso_type == "(1,0)":
            num_inc = count
        elif iso_type == "(0,1)":
            num_exc = count
        elif iso_type == "(1,1)":
            num_both = count

    return num_inc, num_exc, num_both


def filter_event(sample_inc, sample_exc, sample_both,
                 num_total, num_inc, num_exc, num_sum):
    """
    Return True if event passes filter, False otherwise.
    """
    sample_total = sample_inc + sample_exc + sample_both
    sample_sum = sample_inc + sample_exc

    if sample_total < num_total:
        return False

    if sample_sum < num_sum:
        return False

    if sample_inc < num_inc:
        return False

    if sample_exc < num_exc:
        return False

    return True

def multi_filter(filter_filename, output_dir, num_total, num_inc, num_exc,
                 num_sum, delta_psi_filter, bf_filter, vote_thresh,
                 apply_both_samples=False):

    diff_thresh = vote_thresh

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    comp = []
    event_dict = defaultdict(list)
    rep_list = []
    total_events = []

    for fname in filter_filename:

        rep_list.append(defaultdict(list))

        data, h = csv2dictlist_raw(fname)
        total_events.append(len(data))

        filtered = filter_events(data, h, num_total, num_inc, num_exc, num_sum, 
                                    delta_psi_filter, bf_filter, 
                                    apply_both_samples=apply_both_samples)
        comp.append(filtered)

    # if there is just a single file of inputs
    if len(comp) is 1:
        num_pass = len(comp[0])
        output_filename = os.path.join(output_dir,
                                       os.path.basename(fname) + ".filtered")
        print "Filtering %s into %s" %(",".join(filter_filename),
                                       output_filename)
        filter_output(comp[0], output_filename, h, num_pass, total_events[0])

    else:
        # this is not going to work at all, need to figure out how
        # to do this properly
        # create a dictionary to look up the events in
        c = 0
        for events in comp:
            rep_list.append(defaultdict(list))
            for event in events:
                # this dictionary has all the unique event names to check
                event_dict[event['event_name']].append(event)
                # this dictionary has a list of the event names in each replicate
                rep_list[c][event['event_name']].append(event)
            c = c + 1

        for event_name in event_dict.keys():
            # not enough replicates of the event passed the previous filters
            if len(event_dict[event_name]) < diff_thresh:
                del(event_dict[event_name])
                continue

            # test to see if enough of the replicates point in the same direction
            bf_list = []
            dp_list = []
            for event in event_dict[event_name]:
                if not bf_list:
                    bf_pass = bayes_factor_pass(event['bayes_factor'], bf_filter)
                    bf_list = bf_list + bf_pass
                    dp_pass = delta_psi_pass(event['diff'], delta_psi_filter)
                    dp_list = dp_list + dp_pass
                    continue

                # adds the result from the current event to the ones from the 
                # previous events
                bf_pass = bayes_factor_pass(event['bayes_factor'], bf_filter)
                bf_list = [sum(x) for x in zip(bf_list, bf_pass)]
                dp_pass = delta_psi_pass(event['diff'], delta_psi_filter)
                dp_list = [sum(x) for x in zip(dp_list, dp_pass)]

            bf_pass = False
            for bf in bf_list:
                if(abs(bf) >= diff_thresh):
                    bf_pass = True
                    break
            dp_pass = False
            for dp in dp_list:
                if(abs(dp) >= diff_thresh):
                    dp_pass = True
                    break

            if not bf_pass and dp_pass:
                del(event_dict[event_name])

        comp_new = []
        for events in comp:
            event_list = []
            for event in events:
                if event_dict.has_key(event['event_name']):
                    event_list.append(event)
            comp_new.append(event_list)
    
        for i in range(0, len(comp_new)):
            num_pass = len(comp_new[i])
            fname = filter_filename[i]
            output_filename = os.path.join(output_dir,
                                           os.path.basename(fname) + ".filtered")
            print "Filtering %s into %s" %(fname,
                                           output_filename)
            filter_output(comp_new[i], output_filename, h, num_pass, total_events[i])


def bayes_factor_pass(bayes_factor, bf_filter):
    """
    Checks to see out of a list of bayes factors which ones
    pass the filter check. 1 if pass, 0 if no pass.
    """
    if isinstance(bayes_factor, float):
        bayes_factor = [bayes_factor]

    bf_list = []
    for bf in bayes_factor:
        if abs(bf) < bf_filter:
            bf_list.append(0)
        else:
            bf_list.append(1)

    return bf_list

def delta_psi_pass(delta_psi, dp_filter):
    """
    Checks to see out of a list of delta_psi which ones pass the
    filter check. 0 if no pass. -1 if pass but is negative delta psi,
    1 if pass and positive delta psi
    """

    if isinstance(delta_psi, float):
        delta_psi = [delta_psi]

    dp_list = []
    for dp in delta_psi:
        if abs(dp) < dp_filter:
            dp_list.append(0)
        else:
            # this will preserve the sign
            dp_list.append(cmp(dp, 0))
    return dp_list


def fix_bayes_factor(bayes_factor):
    """
    If one of the bayes factors is 'inf' we get a string instead of a 
    tuple back. This is hacky but fixes that.
    """
    # Maximum cut off for Bayes factor value
    max_bf = 1e12
    
    if type(bayes_factor) == str:
        bayes_factor = bayes_factor.split(",")
        bayes_factor = [min(float(x), max_bf) for x in bayes_factor]
        bayes_factor = tuple(bayes_factor)
        bayes_factor = bayes_factor[0]

    return bayes_factor

def filter_events(data, h, num_total, num_inc, num_exc, num_sum,
                  delta_psi_filter, bf_filter,
                  apply_both_samples=False):
    """
    Filter.
    """

    filtered_events = []

    if abs(delta_psi_filter) > 1 or \
       abs(delta_psi_filter) < 0:
        raise Exception, "Error: delta psi value outside [0, 1]." 

    for event in data:
        # Sometimes the bayes factor is not formatted correctly, this fixes that
        event['bayes_factor'] = fix_bayes_factor(event['bayes_factor'])

        num_isoforms = len(event['isoforms'])
        if num_isoforms != 2:
            print "Error: filter_events.py is only defined for MISO output " \
                  "on two-isoform alternative events. " \
                  "Found a non-two isoform event: %s" %(event['event_name'])
            sys.exit(1)

        # Get sample 1 counts
        sample1_counts = get_counts(event['sample1_counts'])

        if sample1_counts == None:
            # Get Bayes factor
            bayes_factor = event['bayes_factor']
            bf_pass = False
            if bayes_factor != 0:
                for bf in bayes_factor:
                    if abs(bf) > abs(bf_filter):
                        bf_pass = True

            if not bf_pass:
                continue

            delta_psi = event['diff']
            dp_pass = False
            for dp in delta_psi:
                if abs(dp) > abs(delta_psi_filter):
                    dp_pass = True

            if not dp_pass:
                continue
        else:
            # Get delta Psi
            delta_psi = float(event['diff'])
            bayes_factor = float(event['bayes_factor'])

            sample1_inc, sample1_exc, sample1_both = sample1_counts

            sample1_result = filter_event(sample1_inc, sample1_exc, sample1_both,
                                          num_total, num_inc, num_exc, num_sum)

            # Get sample 2 counts
            sample2_counts = get_counts(event['sample2_counts'])

            if sample2_counts == None:
                raise Exception, "Incompatible samples."
            
            sample2_inc, sample2_exc, sample2_both = sample2_counts            

            sample2_result = filter_event(sample2_inc, sample2_exc, sample2_both,
                                          num_total, num_inc, num_exc, num_sum)

            if abs(delta_psi) < abs(delta_psi_filter):
                continue

            if abs(bayes_factor) < abs(bf_filter):
                continue

            if apply_both_samples:
                # Apply filter to both samples
                if not sample1_result or not sample2_result:
                    continue
            else:
                if not sample1_result and not sample2_result:
                    continue

        # Event passes filter
        filtered_events.append(event)

    # filter_output(filtered_Events, output_filename, h, num_pass, total_events)
    return(filtered_events)


def filter_output(filtered_events, output_filename, h, num_pass, total_events):
    """
    produce the output of the filtered file
    """
    dictlist2file(filtered_events, output_filename, h)

    print "%d/%d events pass the filter (%.2f percent)." \
          %(num_pass,
            total_events,
            (num_pass/float(total_events)) * 100)


def greeting():
    print "filter_events: filtering MISO pairwise comparison output.\n"
    print "Note: This utility is only works on MISO output for two-isoform "
    print "event annotations.\n"
    
    
def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--filter", dest="filter_filename", default=None, 
                      action="callback", callback=fname_callback,
                      help="Comparison file to filter or list of replicate files to filter.")
    parser.add_option("--control", dest="control_filename", default=[],
                     action="callback", callback=fname_callback,
                     help="Control comparison file to filter.")
    parser.add_option("--output-dir", dest="output_dir", nargs=1,
                      default=None,
                      help="Output directory for filtered file.")
    parser.add_option("--num-total", dest="num_total", nargs=1, default=0,
                      type="int",
                      help="Require at least N-many total reads (inclusion reads, "
                      "exclusion reads or reads supporting both isoforms.")
    parser.add_option("--num-inc", dest="num_inc", nargs=1, default=0,
                      type="int",
                      help="Require at least N-many inclusion reads. "
                      "Assumes that the first isoform (isoform 0) is the "
                      "inclusion isoform.")
    parser.add_option("--num-exc", dest="num_exc", nargs=1, default=0,
                      type="int",
                      help="Require at least N-many exclusion reads")
    parser.add_option("--num-sum-inc-exc", dest="num_sum", nargs=1, default=0,
                      type="int",
                      help="Require the sum of inclusion and exclusion reads to "
                      "be at least N.")
    parser.add_option("--delta-psi", dest="delta_psi", nargs=1, default=0, type="float",
                      help="Require the absolute change in Psi (delta Psi) to be at least N. "
                      "Only takes positive floats as arguments.")
    parser.add_option("--bayes-factor", dest="bayes_factor", nargs=1, default=0, type="float",
                      help="Require the Bayes factor to be at least N.")
    parser.add_option("--apply-both", dest="apply_both", default=False, action="store_true",
                      help="Apply filter to both samples.")
    parser.add_option("--votes", dest="vote_thresh", nargs=1, default=0, type="int",
                    help="The number of biological replicates in a sample which must pass the  "
                    " filters to call an event significant.")

    greeting()

    (options, args) = parser.parse_args()

    if options.filter_filename == None:
        print "Need at least one filename to filter (use --filter.)"
        return
    if options.output_dir == None:
        print "Need an output directory to output filtered file to " \
              "(use --output-dir)"
        return

    filter_filename = []
    for fname in options.filter_filename:
        filter_filename.append(os.path.abspath(os.path.expanduser(fname)))
    control_filename = []
    for fname in options.control_filename:
        control_filename.append(os.path.abspath(os.path.expanduser(fname)))

    output_dir = os.path.abspath(os.path.expanduser(options.output_dir))

    #if len(filter_filename) == 1:
    #    filter_events(filter_filename[1], output_dir,
    #                  options.num_total, options.num_inc, options.num_exc,
    #                  options.num_sum, options.delta_psi, options.bayes_factor,
    #                  apply_both_samples=options.apply_both)

    multi_filter(filter_filename, output_dir, options.num_total, options.num_inc,
            options.num_exc, options.num_sum, options.delta_psi, 
            options.bayes_factor, options.vote_thresh, 
            apply_both_samples=options.apply_both)

if __name__ == '__main__':
    main()
