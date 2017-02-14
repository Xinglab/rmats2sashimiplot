##
## Utilities for working with MISO output samples
##

from scipy import *
from numpy import *

from collections import defaultdict

import os
import sys
import glob
import misopy
import misopy.miso_db as miso_db

from misopy.parse_csv import *
from misopy.credible_intervals import *
import misopy.misc_utils as misc_utils

        
class MISOSamples:
    """
    Representation of MISO samples directory.
    The samples filename is either a plain-text *.miso
    file, e.g. event_name.miso, or a .miso_db file.
    """
    def __init__(self, samples_dir, use_compressed=None):
        self.samples_dir = samples_dir
        # Handle compressed IDS if asked. Load mapping of compressed
        # ids to filenames
        self.compressed_ids_fname = use_compressed
        self.compressed_ids_to_genes = None
        if self.compressed_ids_fname is not None:
            print "  - Loading compressed IDs mapping from: %s" \
                  %(self.compressed_ids_fname)
            # Load mapping from gene IDs to their hashes
            self.compressed_ids_to_genes = \
              misc_utils.load_compressed_ids_to_genes(self.compressed_ids_fname)
            if len(self.compressed_ids_to_genes) == 0:
                print "Error: Compressed IDs shelve file is empty. Are you sure " \
                      "the index directory you passed was created with the " \
                      "--compress-id flag, e.g.:\n" \
                      "index_gff yourfile.gff --compress-id"
                sys.exit(1)
        # Get all the MISO relevant filenames
        self.all_filenames = get_samples_dir_filenames(samples_dir)
        # Mapping of event names to the files that they are in
        # these are either *.miso files or *.miso_db files.
        # Example:
        #   myevent -> myevent.miso
        #   myotherevent_on_chr12 -> chr12.miso_db
        self.event_names_to_fnames = {}
        # Get all the event names in the current samples directory
        self.all_event_names = self.get_all_event_names()
        self.num_events = len(self.all_event_names)

        
    def get_all_event_names(self):
        """
        Return all event names in current samples dir.
        """
        all_event_names = []
        for curr_fname in self.all_filenames:
            if curr_fname.endswith(".miso"):
                # It's a regular .miso plain text file
                event_name = \
                  get_event_name(curr_fname,
                                 use_compressed_map=self.compressed_ids_to_genes)
                # Record event name and its mapping to a .miso file
                all_event_names.append(event_name)
                self.event_names_to_fnames[event_name] = curr_fname
            elif miso_db.is_miso_db_fname(curr_fname):
                # It's a MISO database file, so load all the event
                # names in that file
                curr_db = \
                  miso_db.MISODatabase(curr_fname,
                                       comp_to_uncomp=self.compressed_ids_to_genes)
                # Record event name and its mapping to the chromosome's
                # .miso_db file
                for curr_event_name in curr_db.get_all_event_names():
                    curr_event_name = str(curr_event_name)
                    event_name_to_use = curr_event_name
                    # If we're given a mapping of compressed IDs, use the
                    # mapping to get the uncompressed event name
                    if self.compressed_ids_to_genes is not None:
                        # The internal database representation of compressed
                        # index databases are compressed IDs, so if the
                        # ID is uncompressed it must be converted to a
                        # compressed one.
                        if not misc_utils.is_compressed_name(curr_event_name):
                            event_name_to_use = \
                              str(curr_db.uncomp_to_comp[curr_event_name])
                    all_event_names.append(event_name_to_use)
                    self.event_names_to_fnames[event_name_to_use] = curr_fname
        return all_event_names
    

    def get_event_samples(self, event_name):
        """
        Get the samples information for the given event by name.
        """
        samples = None
        if event_name not in self.event_names_to_fnames:
            return None
        event_fname = self.event_names_to_fnames[event_name]
        if event_fname.endswith(".miso"):
            # Get event from plain text .miso file
            f = open(event_fname, "r")
            samples = load_samples(f)
            f.close()
        elif miso_db.is_miso_db_fname(event_fname):
            # Get event from miso_db file
            curr_db = \
              miso_db.MISODatabase(event_fname,
                                   comp_to_uncomp=self.compressed_ids_to_genes)
            event_data = curr_db.get_event_data_as_stream(event_name)
            samples = load_samples(event_data)
        if samples is None:
            print "WARNING: Could not parse event %s samples" %(event_name)
        return samples
                                

def maxi(l):
    m = max(l)
    for i, v in enumerate(l):
        if m == v:
            return i

        
def load_samples(samples_in):
    """
    Load a file with samples. Takes in a string of data.
    Return the samples, header from the file, the sampled MAP estimate,
    and the sampled MAP's log score.
    """
    try:
        data, h = csv2array(samples_in, skiprows=1, raw_header=True)
        sampled_map_indx = maxi(data['sampled_psi'])
        sampled_map = \
            [float(v) for v in data['sampled_psi'][sampled_map_indx].split(',')]
        sampled_map_log_score = data['log_score'][sampled_map_indx]
        #    print "  - Sampled MAP: %s" %(sampled_map)
        #    print "  - Sampled MAP log score: %.4f" %(sampled_map_log_score)
        samples = []
        for vals in data['sampled_psi']:
            psi_vals = [float(v) for v in vals.split(',')]
            samples.append(psi_vals)
        samples = array(samples)

        # Extract counts from the file's header
        counts_info = get_counts_from_header(h[0])

        return (samples, h, data['log_score'], sampled_map,
                sampled_map_log_score, counts_info)
    except ValueError:
        return None


def parse_sampler_params_from_header(header):
    """
    Parse parameters that were used to produce a set of samples.
    miso_file is a stream to a MISO file.
    """
    if header[0] == '#':
	# strip header start
	header = header[1:]
    fields = header.split('\t')
    params = {}

    for field in fields:
	key, value = field.split('=')
        params[key] = value

    return params


def get_isoforms_from_header(samples_header):
    """
    Given header of a raw MISO samples file, return the isoforms
    field.
    """
    # Get first field (removing prefix comment symbol, '#')
    isoforms = samples_header[1:].split("\t")[0]
    # Remove isoforms= key
    isoforms = isoforms.split("isoforms=")[1]
    # Remove brackets, '[', ']'
    isoforms = isoforms[1:-1]
    
    return isoforms


def get_counts_from_header(samples_header):
    """
    Given a header of a raw MISO samples file, return the
    counts= and assigned_counts= fields.
    """
    fields = samples_header[1:].split("\t")
    counts = {}
    for f in fields:
        if f.startswith("counts="):
            counts['counts'] = f.split("=")[1]
        elif f.startswith("assigned_counts="):
            counts['assigned_counts'] = f.split("=")[1]
            
    if len(counts.keys()) != 2:
        print "Warning: Could not get counts fields out of " \
              "%s header." %(samples_header)
        counts = {'counts': 'n/a',
                  'assigned_counts': 'n/a'}
        
    return counts

    
def get_gene_info_from_params(params):
    """
    Return gene information from parameters of
    a MISO samples file.
    """
    gene_info = defaultdict(lambda: "NA")
    if "chrom" in params:
        gene_info["chrom"] = params["chrom"]
    if "strand" in params:
        gene_info["strand"] = params["strand"]
    if "mRNA_starts" in params:
        gene_info["mRNA_starts"] = params["mRNA_starts"]
    if "mRNA_ends" in params:
        gene_info["mRNA_ends"] = params["mRNA_ends"]
    return gene_info
    

def get_event_name(miso_filename,
                   use_compressed_map=None):
    """
    Get event name from MISO filename.

    Now supports compressed event names.
    """
    basename = os.path.basename(miso_filename)
    if not basename.endswith(".miso"):
        # Not a MISO filename
        return None
    event_name = basename.split(".miso")[0]
    if use_compressed_map is not None:
        if event_name not in use_compressed_map:
            print "MISO FILENAME IS: %s" %(miso_filename)
            print event_name
        else:
            event_name = use_compressed_map[event_name]
    return event_name

# def get_event_name(miso_filename):
#     """
#     Get event name from MISO filename.
#     """
#     basename = os.path.basename(miso_filename)
#     if not basename.endswith(".miso"):
#         # Not a MISO filename
#         return None
#     event_name = basename.split(".miso")[0]
#     return event_name
    
    
def summarize_sampler_results(samples_dir, summary_filename,
                              use_compressed=None):
    """
    Given a set of samples from MISO, output a summary file.
    """
    summary_file = open(summary_filename, 'w')
    header_fields = ["event_name", "miso_posterior_mean", "ci_low", "ci_high",
                     "isoforms", "counts", "assigned_counts",
                     # Fields related to gene/event
                     "chrom",
                     "strand",
                     "mRNA_starts",
                     "mRNA_ends"]
    summary_header = "%s\n" %("\t".join(header_fields))
    summary_file.write(summary_header)
    print "Loading events from: %s" %(samples_dir)
    print "Writing summary to: %s" %(summary_filename)
    samples_obj = MISOSamples(samples_dir,
                              use_compressed=use_compressed)
    num_events = 0

    for event_name in samples_obj.all_event_names:
        samples_results = samples_obj.get_event_samples(event_name)
        if samples_results is None:
            print "WARNING: Skipping %s" %(event_name)
            # Skip files that could not be parsed
            continue
        # If we're not given a mapping to compressed IDs, check
        # that the event IDs do not look compressed
        if misc_utils.is_compressed_name(event_name) and \
           (use_compressed is None):
            print "WARNING: %s looks like a compressed id, but no mapping file " \
                  "from compressed IDs to event IDs was given! Try: --use-compressed" \
                  %(event_name)
        # Load header/parameters information
        samples = samples_results[0]
        header = samples_results[1]
        header = header[0]
        params = parse_sampler_params_from_header(header)
        # Get counts information from header
        counts_info = samples_results[5]
        shape_len = len(shape(samples))
        if shape_len < 2:
            print "WARNING: Skipping %s -- mishaped file" %(event_name)
            continue
        num_samples, num_isoforms = shape(samples)
        output_fields = format_credible_intervals(event_name, samples)
            
        # Add isoforms information to output fields
        isoforms_field = get_isoforms_from_header(header)
        output_fields.append(isoforms_field)

        # Add counts information to output fields
        output_fields.append(counts_info['counts'])
        output_fields.append(counts_info['assigned_counts'])

        gene_info = get_gene_info_from_params(params)
        output_fields.append(gene_info["chrom"])
        output_fields.append(gene_info["strand"])
        output_fields.append(gene_info["mRNA_starts"])
        output_fields.append(gene_info["mRNA_ends"])
        
        output_line = "%s\n" %("\t".join(output_fields))
	summary_file.write(output_line)
	num_events += 1
    print "  - Summarized a total of %d events." %(num_events)
    summary_file.close()

    
def is_miso_chrom_dir(dirname):
    """
    Return True if a directory contains *.miso files.
    """
    if not os.path.isdir(dirname):
        return False
    basename = os.path.basename(dirname)
    # If the directory is named like a chromosome directory,
    # keep it
    if basename.startswith("chr") or basename.isdigit() or \
        basename == "X" or basename == "Y":
        return True
    # If naming is unclear, check that it contains *.miso files
    fnames = glob.glob(os.path.join(dirname, "*.miso"))
    if len(fnames) >= 1:
        return True
    return False
    
    
def get_samples_dir_filenames(samples_dir):
    """
    Get all the filenames associated with a samples directory.

    Assumes samples directory have the following structure:

      - samples_dir
        - chr1
        - chr2
        ...
        - chrN

    or:

      - samples_dir
        - chr1.miso_db
        - chr2.miso_db
        ...
        - chrN.miso_db

    Also collect files in samples_dir for backwards compatibility.
    """
    directories = glob.glob(os.path.join(samples_dir, "*"))
    directories = filter(is_miso_chrom_dir, directories)
    
    # Filenames indexed by chromosomes
    filenames = []

    for directory in directories:
        if os.path.isdir(directory):
            dir_filenames = os.listdir(directory)
            dir_filenames = [os.path.join(directory, dname) \
                             for dname in dir_filenames]
            filenames.extend(dir_filenames)

    # Filenames in top-level directory
    filenames.extend(os.listdir(samples_dir))

    # Add parent directory to all filenames
    filenames = [os.path.join(samples_dir, fname) \
                 for fname in filenames]

    # Remove directories and files beginning with "."
    filenames = filter(lambda f: not os.path.isdir(f),
                       filenames)
    filenames = filter(lambda f: not os.path.basename(f).startswith("."),
                       filenames)

    # Resulting files should be either *.miso files
    # or *.miso_db files, but not both
    miso_filenames = \
      filter(lambda f: os.path.basename(f).endswith(".miso"),
             filenames)
    miso_db_filenames = \
      filter(miso_db.is_miso_db_fname,
             filenames)
    if len(miso_filenames) > 0 and len(miso_db_filenames) > 0:
        print "WARNING: Directory %s has both *.miso and *.miso_db files" \
              %(samples_dir)
    relevant_filenames = miso_filenames + miso_db_filenames
    return relevant_filenames


