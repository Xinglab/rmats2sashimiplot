##
## Representation for AS events
##
import misopy
import misopy.json_utils as json_utils
import misopy.pickle_utils as pickle_utils
import misopy.Gene as Gene
import time
import csv
import os

class Event(object):
    """
    Representation of an alternative splicing event.
    """
    def __init__(self, label, pcr_psi=None, uniqueness=None,
                 chrom=None):
	self.label = label
	# The gene that it's in
	#self.gene = gene
	self.pcr_psi = pcr_psi
	self.psi_posterior_mean = None
	# representation of uniqueness
	self.uniqueness = uniqueness
        self.chrom = chrom

class TwoIsoEvent(Event):
    """
    Representation of an alternative splicing event with two isoforms.
    """
    def __init__(self, label, event_type, chrom,
		 psi_map=None,
		 psi_sj=None,
		 pcr_psi=None,
		 psi_posterior_mean=None,
		 uniquness=None):
	Event.__init__(self, label, chrom=chrom)
	# Event's length information
	# length of AS part
	self.label = label
	self.event_type = event_type
	self.len = None
	self.up_part_len = None
	self.dn_part_len = None
	
	# Representation of RNA-Seq counts
	# for SE events
	self.num_inc = None
	self.num_exc = None
	self.num_common = None
	# for TandemUTR events
	self.num_core = None
	self.num_ext = None
        # for AFE/ALE events
        self.num_proximal_body = None
        self.num_distal_body = None
        self.num_proximal_jxns = None
        self.num_distal_jxns = None
	
	# Estimates of Psi values
	self.psi_map = psi_map
	self.psi_sj = psi_sj
	self.pcr_psi = pcr_psi
	self.psi_posterior_mean = psi_posterior_mean

    def __repr__(self):
	if self.event_type == 'SE':
	    return "TwoIsoEvent(%s, %s)[ni = %d, ne = %d, nb = %d]" %(self.label,
								      self.event_type,
								      self.num_inc,
								      self.num_exc,
								      self.num_common)
	elif self.event_type == 'TandemUTR':
	    return "TwoIsoEvent(%s, %s)[num_core = %d, num_ext = %d]" %(self.label,
									self.event_type,
									self.num_core,
									self.num_ext)
        elif self.event_type == 'RI':
	    return "TwoIsoEvent(%s, %s)[ni = %d, ne = %d, nb = %d]" %(self.label,
                                                                      self.event_type,
                                                                      self.num_inc,
                                                                      self.num_exc,
                                                                      self.num_common)
	else:
	    return "UnsupportedTwoIsoEvent(%s, %s)" %(self.label,
                                                      self.event_type)

class MultIsoEvent(Event):
    """
    Representation of an alternative splicing event with multiple isoforms.
    """
    def __init__(self, label, gene):
	Event.__init__(self, label)
	self.gene = gene
	
    def __repr__(self):
	return "MultIsoEvent(label = %s)" %(self.label)

class MISOEvents:
    """
    A set of AS events for use by MISO.
    """
    def __init__(self, num_iso, event_type, events=None, from_file=None):
	self.num_iso = num_iso
	self.event_type = event_type
	# Can't give a set of events and load events from a file at the same time
	assert(not(events and from_file))
	self.events = events
        self.from_file = from_file
	if self.from_file != None:
            self.from_file = os.path.expanduser(from_file)            
	    #self.load_from_json_file(self.from_file)
	    self.load_from_pickle_file(self.from_file)
	self.num_events = len(self.events)


#     def get_event(self, event_name):
# 	for event in self.events:
# 	    if event.label == event_name:
# 		return event
# 	return None

    def get_event(self, event_name):
	"""
	Assume each event has a unique identifier.
	"""
	return self.events[event_name]

    def filter_events(self, settings=None):
	"""
	Remove events that do not meet coverage criteria.
	"""
	if self.event_type == 'SE':
	    self.filter_se_events(settings)
        elif self.event_type == 'RI':
            self.filter_ri_events(settings)
	elif self.event_type == 'TandemUTR':
	    self.filter_tandem_utr_events(settings)
	elif (self.event_type == 'AFE' or self.event_type == 'ALE'):
            self.filter_afe_ale_events(settings)
        else:
	    raise Exception, "Unsupported event type for filtering: %s" %(self.event_type)

    def filter_afe_ale_events(self, settings,
                              atleast_proximal=0,
                              atleast_distal=0,
                              proximal_distal_sum=20):
        print "Filtering AFE/ALE events..."
        filtered_events = {}
        for event_name, event in self.events.iteritems():
            num_proximal = event.num_proximal_body + event.num_proximal_jxns
            num_distal = event.num_distal_body + event.num_distal_jxns
            if (num_proximal >= atleast_proximal and num_distal >= atleast_distal) \
                and (num_proximal + num_distal) >= proximal_distal_sum:
                filtered_events[event_name] = event
        self.events = filtered_events

    def filter_tandem_utr_events(self, settings,
                                 atleast_core=1,
                                 atleast_ext=1,
                                 core_ext_sum=20):
	print "Filtering tandem UTR events..."
        if settings != None:
            if 'utr_filter' in settings:
                core_ext_sum = settings['utr_filter'][0]
                atleast_ext = settings['utr_filter'][1]
                atleast_core = settings['utr_filter'][2]
        
	filtered_events = {}
	for event_name, event in self.events.iteritems():
	    if (event.num_core >= atleast_core and event.num_ext >= atleast_ext) and \
	       (event.num_core + event.num_ext) >= core_ext_sum:
		filtered_events[event_name] = event
	self.events = filtered_events

    def filter_ri_events(self, settings,
                         atleast_ri_plus_ne=10,
                         atleast_ne=0,
                         atleast_num_common=1):
        print "Filtering RI events..."
        print "Filter: "
        print "  - ri_plus_ne >= %d" %(atleast_ri_plus_ne)
        print "  - ne >= %d" %(atleast_ne)
        print "  - num_common >= %d" %(atleast_num_common)
            
	filtered_events = {}
        for event_name, event in self.events.iteritems():
	    if ((event.num_inc + event.num_exc) >= atleast_ri_plus_ne \
                and (event.num_exc >= atleast_ne) and \
                (event.num_common >= atleast_num_common)):
		filtered_events[event_name] = event
	self.events = filtered_events

    def filter_se_events(self, settings,
                         atleast_ni_plus_ne=10,
                         atleast_ne=0,
                         atleast_num_common=1):
	print "Filtering SE events..."
        if settings != None:
            if 'se_filter' in settings:
                atleast_ni_plus_ne = settings['se_filter'][0]
                atleast_ne = settings['se_filter'][1]
                atleast_num_common = settings['se_filter'][2]
            
	filtered_events = {}
        for event_name, event in self.events.iteritems():
	    if ((event.num_inc + event.num_exc) >= atleast_ni_plus_ne) and (event.num_exc >= atleast_ne) and \
		   (event.num_common >= atleast_num_common):
		filtered_events[event_name] = event
	self.events = filtered_events

    def load_from_json_file(self, json_filename):
	# clear currently loaded events, if any
	self.clear_events()
	# Modify events directly
	t1 = time.time()
	self.events = json_utils.json_load_file(json_filename)
	t2 = time.time()
	print "Loading from JSON file took %.2f seconds." %(float(t2 - t1))
	self.num_events = len(self.events)

    def load_from_pickle_file(self, pickle_filename):
        print "Called on: ", pickle_filename
	# clear currently loaded events, if any
	self.clear_events()
	# Modify events directly
	t1 = time.time()
	self.events = pickle_utils.load_pickled_file(pickle_filename)
	t2 = time.time()
	print "Loading from Pickle file took %.2f seconds." %(float(t2 - t1))
	self.num_events = len(self.events)

    def loaded_events_to_genes(self, single_event_name=None,
                               read_len=None, overhang_len=None):
	"""
	Parse the loaded set of events into gene structures.  Map events to genes.
	"""
	if len(self.events) == 0:
	    raise Exception, "Must load events first before they can be converted to genes."
	events_to_genes = {}

	t1 = time.time()
	if single_event_name:
	    # If given an event name, only parse that event
	    event_names = [single_event_name]
	else:
	    event_names = self.events.keys()
	for event_name in event_names:
	    event = self.events[event_name]

	    if self.event_type == 'SE' or self.event_type == 'RI':
		gene = Gene.se_event_to_gene(event.up_part_len, event.len, event.dn_part_len,
                                             event.chrom,
                                             label=event.label)
	    elif self.event_type == 'TandemUTR':
		gene = Gene.tandem_utr_event_to_gene(event.core_len, event.ext_len,
                                                     event.chrom,
                                                     label=event.label)
            elif (self.event_type == 'AFE' or self.event_type == 'ALE'):
                gene = Gene.afe_ale_event_to_gene(event.proximal_exons, event.distal_exons,
                                                  self.event_type, event.chrom,
                                                  label=event.label, read_len=read_len,
                                                  overhang_len=overhang_len)
            else:
                raise Exception, "Unsupported event type: %s" %(self.event_type)
	    events_to_genes[event_name] = gene
	t2 = time.time()
	print "Parsing of events to genes took %.2f seconds." %(t2 - t1)
	return events_to_genes

    def output_file(self, results_output_dir, sample_label, method="pickle"):
	events_filename = None
	if method == "json":
	    events_filename = self.output_json_file(results_output_dir, sample_label)
	elif method == "pickle":
	    events_filename = self.output_pickle_file(results_output_dir, sample_label)
	return events_filename

    def output_json_file(self, results_output_dir, sample_label):
	"""
	Output as json.
	"""
	print "Serializing a total of %d events by JSON." %(len(self.events))
	json_output_dir = os.path.join(results_output_dir, 'json')
	if not os.path.isdir(json_output_dir):
	    os.mkdir(json_output_dir)
	json_events_filename = os.path.join(json_output_dir, sample_label + '.json')
	json_utils.json_serialize(self.events, json_events_filename)
	return json_events_filename

    def output_pickle_file(self, results_output_dir, sample_label):
	print "Serializing a total of %d events by Pickle." %(len(self.events))
	pickle_output_dir = os.path.join(results_output_dir, 'pickle')
	if not os.path.isdir(pickle_output_dir):
	    os.mkdir(pickle_output_dir)
	pickle_events_filename = os.path.join(pickle_output_dir, sample_label + '.pickle')
	pickle_utils.write_pickled_file(self.events, pickle_events_filename)
	return pickle_events_filename

    def clear_events(self):
	self.events = []
    
def parse_part(exon, delimiter=':'):
    chrom, start_coord, end_coord, strand = exon.split(delimiter)
    start_coord = int(start_coord)
    end_coord = int(end_coord)
    exon_len = abs(end_coord - start_coord) + 1
    return {'chrom': chrom,
	    'start_coord': start_coord,
	    'end_coord': end_coord,
	    'strand': strand,
	    'len': exon_len}

def parse_event_information(event_name, event_type, delimiter=';',
                            events_to_info=None):
    """
    Parse event information from event counts file.
    """
    if event_type == 'SE' or event_type == 'RI':
	# Parse lengths of exons in event
	up_part, se_part, dn_part = event_name.split(delimiter)
	up_part_info = parse_part(up_part)
	se_part_info = parse_part(se_part)
	dn_part_info = parse_part(dn_part)
	return {'up_part': up_part_info,
		'part': se_part_info,
		'dn_part': dn_part_info,
                'chrom': up_part_info['chrom']}
    elif event_type == 'TandemUTR':
	core_part, ext_part = event_name.split(delimiter)
	core_part_info = parse_part(core_part)
	ext_part_info = parse_part(ext_part)
	return {'core_part': core_part_info,
		'ext_part': ext_part_info,
                'chrom': core_part_info['chrom']}
    elif (event_type == 'AFE' or event_type == 'ALE'):
        if event_name not in events_to_info:
            raise Exception, "Error: Given unknown event %s of type %s." \
                  %(event_name, event_type)
        # Return information about the events
        return events_to_info[event_name]

def parse_afe_ale_event(proximal_exons_str, distal_exons_str,
                        delimiter=','):
    """
    Extract lengths of AFE/ALE proximal and distal exons.
    """
    proximal_exons = proximal_exons_str.split(delimiter)
    distal_exons = distal_exons_str.split(delimiter)

    event_info = {}

    proximal_exons = [parse_part(proximal_exon) \
                      for proximal_exon in proximal_exons]
    distal_exons = [parse_part(distal_exon) \
                    for distal_exon in distal_exons]

    assert(len(proximal_exons) > 0)
    assert(len(distal_exons) > 0)
    
    event_info = {'proximal_exons':
                  proximal_exons,
                  'distal_exons':
                  distal_exons}
    return event_info
    
def load_afe_ale_events_information(events_info_filename, event_type,
                                    delimiter='\t'):
    """
    Load information about AFE/ALE events.
    """
    assert((event_type == 'AFE') or (event_type == 'ALE')), \
                       "Error: Event type must be AFE/ALE"

    print "Loading events from %s (event type: %s)" %(events_info_filename,
                                                      event_type)
    events_info_file = open(events_info_filename, 'r')
    events_info = csv.reader(events_info_file,
                             delimiter=delimiter)
    events_to_info = {}
    
    for event_name, proximal_exons, distal_exons in events_info:
        events_to_info[event_name] = parse_afe_ale_event(proximal_exons,
                                                         distal_exons)
    events_info_file.close()
    return events_to_info
    
def load_event_counts(events_filename, event_type, delimiter=';',
		      events_info_filename=None):
    """
    Parse mRNA-Seq counts data.
    """
    events_file = open(events_filename)
    events = {}
    events_to_info = None

    # Load length information about events, if given
    if events_info_filename != None:
        if (event_type == 'AFE' or event_type == 'AFE'):
            events_to_info = load_afe_ale_events_information(events_info_filename,
                                                             event_type)
        
    for line in events_file:
	event_name, counts = line.strip().split('\t')
	counts_list = counts.split(delimiter)
	
	# Make sure there's more than one count in the set
	assert(len(counts_list) > 1)

	counts = [int(c) for c in counts_list]

        # Part information about event (lengths of exons, etc.)
	event_info = parse_event_information(event_name, event_type,
                                             events_to_info=events_to_info)

        chrom = event_info['chrom']

	event = None

        ##
        ## Skipped exons
        ##
	if event_type == 'SE':
	    up_part_len = event_info['up_part']['len']
	    part_len = event_info['part']['len']
	    dn_part_len = event_info['dn_part']['len']

	    # Parse counts
	    num_up, num_se, num_dn, num_upinc, num_dninc, num_exc = counts
	    num_inc = num_se + num_upinc + num_dninc
	    num_common = num_up + num_dn

	    # Compile into event
	    event = TwoIsoEvent(event_name, event_type, chrom)
	    event.len = part_len
	    event.up_part_len = up_part_len
	    event.dn_part_len = dn_part_len
	    event.num_inc = num_inc
	    event.num_exc = num_exc
	    event.num_common = num_common
        ##
        ## Tandem UTR
        ##
	elif event_type == 'TandemUTR':
	    # Parse counts
            # Assume the format is: core; ext
	    num_ext, num_core = counts

	    # Compile into event
	    event = TwoIsoEvent(event_name, event_type, chrom)
	    event.core_len = event_info['core_part']['len']
	    event.ext_len = event_info['ext_part']['len']
	    
	    event.num_core = num_core
	    event.num_ext = num_ext
        ##
        ## AFE/ALE
        ##
	elif (event_type == 'AFE' or event_type == 'ALE'):
            # Format of counts for AFE/ALE is:
            # counts for: proximal body, distal body, proximal junctions, and distal junctions
            num_proximal_body, num_distal_body, num_proximal_jxns, num_distal_jxns = counts
            
            event = TwoIsoEvent(event_name, event_type, chrom)
            event.proximal_exons = event_info['proximal_exons']
            event.distal_exons = event_info['distal_exons']
            event.num_proximal_body = num_proximal_body
            event.num_distal_body = num_distal_body
            event.num_proximal_jxns = num_proximal_jxns
            event.num_distal_jxns = num_distal_jxns
        ##
        ## Retained intron
        ##
	elif event_type == 'RI':
            # Format of RI counts
            # up;ri;dn;ejxn
            num_up, num_ri, num_dn, num_exc = counts
            event = TwoIsoEvent(event_name, event_type, chrom)
            event.num_inc = num_ri
            event.num_exc = num_exc
            event.num_common = num_up + num_dn
            
            # Length information
	    up_part_len = event_info['up_part']['len']
	    part_len = event_info['part']['len']
	    dn_part_len = event_info['dn_part']['len']

	    event.len = part_len
	    event.up_part_len = up_part_len
	    event.dn_part_len = dn_part_len
        ##
        ## MXEs
        ##
	elif event_type == 'MXE':
	    raise Exception, "MXEs not supported."

	assert(event != None), "Event type %s is unknown." %(event_type)
	    
	events[event_name] = event
	
    miso_events = MISOEvents(2, event_type, events=events)
    return miso_events
