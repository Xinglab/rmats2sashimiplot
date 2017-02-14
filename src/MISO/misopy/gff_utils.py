##
## Revamped GFF utils
##
# Copyright (c) 2007
# Colin Dewey (University of Wisconsin-Madison)
# cdewey@biostat.wisc.edu
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# References for GFF formats
# version   reference
# 1 & 2     http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml
# 2.1       http://genes.cs.wustl.edu/GTF21.html
# 2.2       http://genes.cs.wustl.edu/GTF22.html
# 2.5       ???
# 3         http://song.sourceforge.net/gff3.shtml

import os
import sys
import re
import shelve
import misopy
import misopy.pickle_utils as pickle_utils
from urllib import quote as url_quote, unquote as url_unquote

from collections import defaultdict

#__all__ = ["GFF", "GFFDatabase", "Reader", "Writer", "FormatError",
#           "Metadatum", "SequenceRegion"]

# def get_gene_record_ids(gff_filename, version="3"):
#     """
#     Pull all the gene IDs from a GFF.
#     """
#     FILE = open(gff_filename, "r")
#     gff = Reader(FILE, version)

#     gene_records = [record for record in gff \
#                     if record.type == "gene"]

#     return gene_records

def load_indexed_gff_file(indexed_gff_filename):
    """
    Load indexed representation of a set of genes.
    """
    indexed_gff = pickle_utils.load_pickled_file(indexed_gff_filename)
    return indexed_gff

    
def load_indexed_gff_chrom(indexed_gff_chrom_filename):
    """
    Load indexed representation of a GFF chromosome.
    """
    indexed_gff_chrom = pickle_utils.load_pickled_file(indexed_gff_chrom_filename)
    return indexed_gff_chrom


def load_shelved_genes_to_fnames(indexed_gff_dir,
                                 shelve_basename="genes_to_filenames.shelve"):
    """
    Load mapping from gene IDs to their indexed
    filenames from the 'genes_to_filenames.shelve' file
    if it exists. Return None if it does not exist.
    """
    shelve_fname = os.path.join(indexed_gff_dir, shelve_basename)
    print "Searching for %s.." %(shelve_fname)
    gene_ids_to_gff_index = None    
    if os.path.isfile(shelve_fname):
        print "  - Found shelved file."
        gene_ids_to_gff_index = shelve.open(shelve_fname)
    else:
        print "  - File not found."
    return gene_ids_to_gff_index


def get_gene_ids_to_gff_index(indexed_gff_dir):
    """
    Return mapping from gene IDs to their indexed GFF filename.
    """
    print "Mapping genes to their indexed GFF representation, using %s" \
          %(indexed_gff_dir)
    gff_chrom_dirs = os.listdir(indexed_gff_dir)

    # Load gene IDs to mapping from .shelve file
    # if it exists
    gene_ids_to_gff_index = load_shelved_genes_to_fnames(indexed_gff_dir)
    if gene_ids_to_gff_index is not None:
        # File loaded from shelve successfully
        return gene_ids_to_gff_index

    # Shelve file not found, so reconstruct the mapping from gene IDs
    # to their indexed filenames
    gene_ids_to_gff_index = {}

    for chrom_dir in gff_chrom_dirs:
        chrom_dir_path = os.path.abspath(os.path.join(indexed_gff_dir,
                                                      chrom_dir))

        # Skip subentries that are not directories
        if not os.path.isdir(chrom_dir_path):
            print "Skipping: %s" %(chrom_dir_path)
            continue

        # Get the chromosome filename
        chrom_indexed_filenames = os.listdir(chrom_dir_path)

        # Check if the indexed representation is one file
        # corresponding to a chromosome or multiple files
        # corresponding to many genes
        num_genes = len(chrom_indexed_filenames)
        if num_genes > 1:
            # Handle genes case
            print "Loading indexed gene filenames from: %s" \
                  %(chrom_dir_path)
            print "  - Loading %d genes" %(num_genes)
            
            for gene_index_filename in chrom_indexed_filenames:
                # Skip non-Pickle files
                if not gene_index_filename.endswith(".pickle"):
                    continue

                gene_index_filename = os.path.abspath(os.path.join(chrom_dir_path,
                                                                   gene_index_filename))
                indexed_gene = load_indexed_gff_chrom(gene_index_filename)
                
                for gene_id, gene_info in indexed_gene.iteritems():
                    gene_ids_to_gff_index[gene_id] = gene_index_filename
        elif num_genes == 0:
            raise Exception, "No genes in directory: %s" %(chrom_dir_path)
        else:
            chrom_indexed_filename = os.path.join(chrom_dir_path,
                                                  chrom_indexed_filenames[0])

            # Load the indexed chromosome GFF file
            indexed_gff_chrom = load_indexed_gff_chrom(chrom_indexed_filename)

            for gene_id, gene_info in indexed_gff_chrom.iteritems():
                gene_ids_to_gff_index[gene_id] = chrom_indexed_filename

    return gene_ids_to_gff_index


def parse_gff_attribs(attrib_str):
    attribs = {}
    for pair in attrib_str.split(";"):
        key, val = pair.split("=")
        attribs[key] = val
    return attribs
                
        
class GFFDatabase:
    """
    A set of GFF entries from a GFF file.
    """
    def __init__(self, from_filename=None,
                 reverse_recs=False,
                 include_introns=False,
                 suppress_warnings=False):
	self.genes = []
	self.mRNAs = []
	self.exons = []
	self.cdss = []
        self.__entries = []
	self.from_filename = from_filename
        self.suppress_warnings = suppress_warnings

        ## Indexed representation of GFFs
        self.mRNAs_by_gene = defaultdict(list)
        self.exons_by_mRNA = defaultdict(list)
        self.cdss_by_exon = defaultdict(list)

        self.suppress_warnings = suppress_warnings

	if from_filename:
	    # load GFF from given filename
	    self.from_file(from_filename,
                           reverse_recs=reverse_recs,
                           include_introns=include_introns)
	    self.from_filename = from_filename
            
    def __len(self):
        return len(self.__entries)


    def from_file(self, filename, version="3",
                  reverse_recs=False,
                  include_introns=False):
        FILE = open(filename, "r")
        reader = Reader(FILE, version)
        for record in reader.read_recs(reverse_recs=reverse_recs):
	    if record.type == "gene":
		self.genes.append(record)
            # Allow "transcript" 
	    elif record.type == "mRNA" or record.type == "transcript":
		self.mRNAs.append(record)
                self.mRNAs_by_gene[record.get_parent()].append(record)
	    elif record.type == "exon":
		self.exons.append(record)
                self.exons_by_mRNA[record.get_parent()].append(record)
            elif include_introns and (record.type == "intron"):
                # Treat introns like exons if asked to
                self.exons.append(record)
                self.exons_by_mRNA[record.get_parent()].append(record)
	    elif record.type == "CDS":
		self.cdss.append(record)
                self.cdss_by_exon[record.get_parent()].append(record)
	    # is there a need to store all entries separately? Probably not but
	    # leaving it in for now
            self.__entries.append(record)
	self.from_filename = filename
        FILE.close()

    def get_genes_records(self, genes):
	"""
	Return all the relevant records for a set of genes.
	"""
	recs = []
        gene_hierarchy = {}
        
	for gene in genes:
	    mRNAs = []
	    exons = []
	    cdss = []

            # Initialize hierarchical structure per gene
            gene_hierarchy[gene] = {'mRNAs': defaultdict(dict)}

            genes_mRNAs = self.mRNAs_by_gene[gene]
            
	    # find all the relevant mRNAs
            for mRNA_rec in genes_mRNAs:
                mRNA_rec_id = mRNA_rec.get_id()
                # Initialize structure per mRNA
                gene_hierarchy[gene]['mRNAs'][mRNA_rec_id] = \
                                                           {'exons': defaultdict(dict),
                                                            'record': mRNA_rec}
                mRNAs.append(mRNA_rec)

            # Find all the gene's exons
            all_exons_of_gene = [self.exons_by_mRNA[mrna] \
                                 for mrna in genes_mRNAs]
            
	    # find all the exons of each of the gene's mRNAs
            for mRNA_rec in genes_mRNAs:
                mRNA_rec_id = mRNA_rec.get_id()
                genes_exons = self.exons_by_mRNA[mRNA_rec_id]
                
                for exon_rec in genes_exons:
                    mRNA_rec_id = mRNA_rec.get_id()
                    exon_rec_id = exon_rec.get_id()
                    gene_hierarchy[gene]['mRNAs'][mRNA_rec_id]['exons'][exon_rec_id] = \
                                  {'cdss': defaultdict(list),
                                   'record': exon_rec}
                    exons.append(exon_rec)

                    # for each exon, find the cdss
                    exon_cdss = self.cdss_by_exon[exon_rec_id]
                    for cds_rec in exon_cdss:
                        cds_rec_id = cds_rec.get_id()                    
                        gene_hierarchy[gene]['mRNAs'][mRNA_rec_id]['exons'][exon_rec_id]['cdss'][cds_rec_id] = \
                                                      {'record': cds_rec}
			cdss.append(cds_rec)
                        
	    if len(mRNAs) == len(exons) == len(cdss) == 0:
                if not self.suppress_warnings:
                    print "WARNING: No entries found for gene %s in GFF %s" \
                          %(gene, self.from_filename)
                # Remove from gene hierarchy
                del gene_hierarchy[gene]
#		raise Exception, "No entries found for gene %s in GFF: %s" %(gene,
#                                                                             self.from_filename)
	    # add mRNAs
	    recs.extend(mRNAs)
	    # add exons
	    recs.extend(exons)
	    # add cdss
	    recs.extend(cdss)
	return recs, gene_hierarchy
    

    def write_genes(self, genes, filename):
	"""
	Serialize a set of genes (list of IDs for the genes field) to the given filename as GFF.
	"""
	recs, hierarchy = self.get_genes_records(genes)
	# retrieve the records corresponding to the genes
	if len(recs) == 0:
	    raise Exception, "No entries found for " + str(genes) + " in GFF: %s" %(filename)
	# serialize them as GFF
	output_file = open(filename, 'w')
	print >> sys.stderr, "Outputting sliced GFF records to: %s" %(filename)
	gff_writer = Writer(output_file)
	gff_writer.write_recs(recs)

    def next(self):
	if self.__entries == []:
	    raise StopIteration
        return self.__entries.pop()

    def __iter__(self):
	return self

class GFF:
    """A record from a GFF file.

    Fields:
       seqid
       source
       type
       start
       end
       score
       strand
       phase
       attributes
    """

    def __init__(self, seqid, source, type, start, end,
                 score=None, strand=None, phase=None, attributes=None):
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase

        if attributes:
            self.attributes = attributes
        else:
            self.attributes = {}

        # sanitize: require that start <= end
        if self.start > self.end:
            self.start, self.end = self.end, self.start
            if strand != '-':
                print >>sys.stderr, "WARNING: Swapping start and end fields, which must satisfy start <= end:"
                print >>sys.stderr, self.__repr__()

        # set default exon IDs if they are not already set
        self._set_default_exon_id()

        # Filter the IDs to not have underscores, since these are used
        # internally.
        self._filter_exon_id()
        

    def _set_default_exon_id(self):
        """
        If exon records are missing an ID, automatically compute an ID of the form:

        parent_transcript@start_coord@end_coord@strand

        This function is dedicated to Robert K. Bradley.
        """
        if self.type == "exon" and "ID" not in self.attributes:
            parent_id = self.get_parent()
            exon_id = "%s@%s@%s@%s" %(parent_id,
                                      self.start,
                                      self.end,
                                      self.strand)
            self.attributes['ID'] = [exon_id]

    def _filter_exon_id(self, replace_char='@'):
        """
        Replace exon ID underscores (_) with another symbol.
        """
        return
        # delim_symbol = "_"
        # if self.type == "exon" and "ID" in self.attributes:
        #     exon_id = self.get_id()
        #     if delim_symbol in exon_id:
        #         new_exon_id = exon_id.replace(delim_symbol, replace_char)
        #         print "Warning: replacing exon id %s with %s" %(exon_id,
        #                                                         new_exon_id)
        #         self.attributes['ID'] = [new_exon_id]
            

    def copy(self):
        """Returns a copy of this GFF record"""

        # Make new attributes dict with copied value lists
        attributes_copy = dict([(k, v[:]) for k, v in self.attributes.items()])

        return GFF(self.seqid,
                      self.source,
                      self.type,
                      self.start,
                      self.end,
                      score=self.score,
                      strand=self.strand,
                      phase=self.phase,
                      attributes=attributes_copy)

    def is_valid(self):
        """Returns True if this record passes basic GFF record requirements."""

        # Check that the first five fields are defined and non-empty/zero
        if not (self.seqid and self.source and self.type and
                self.start and self.end):
            return False

        # Check that [start, end] is a valid genomic interval
        if not (is_integer(self.start) and is_integer(self.end) and
                  self.start > 0 and self.start <= self.end):
            return False
        
        # Check that score is a valid floating point number
        if self.score is not None:
            try:
                float(self.score)
            except ValueError:
                return False
        
        if self.strand not in (None, '+', '-'):
            return False
        if self.phase not in (None, 0, 1, 2):
            return False
        
        # Check that CDS records have phase defined
        if self.type == "CDS" and self.phase is None:
            return False

        return True

    def __repr__(self):
        return "GFF(%s, %s, %s, %s, %s, %s, %s, %s, %s)" % \
            tuple(map(repr, (self.seqid, self.source, self.type, 
                             self.start, self.end,
                             self.score, self.strand, self.phase, 
                             self.attributes)))

    def __len(self):
        """Compute the length of the feature."""

        return self.end - self.start + 1

    def length(self):
        """Get the length of a feature."""
        return self.__len()

    def get_values(self, key):
        """Get the values of a particular key.  Return the empty list if no such key."""

        if key in self.attributes:
            return self.attributes[key]
        return []

    def get_value(self, key):
        """Get the value of a particular key.  If multiple values, return the first value.
        If no such key, return the empty list."""

        if key in self.attributes:
            # Ensure that trailing whitespace is removed
            return self.attributes[key][0].rstrip()
        return ""
 
    def get_id (self):
        """Get the ID attribute."""
        return self.get_value ("ID")

    def get_parent (self):
        """Get the Parent attribute."""
        
        return self.get_value ("Parent")

    def get_name (self):
        """Get the Name attribute."""
        
        return self.get_value ("Name")

    def get_note (self):
        """Get the Note attribute."""
        
        return self.get_value ("Note")

class Metadatum:
    def __init__(self, name, value=None):
        self.name = name
        self.value = value

class SequenceRegion(Metadatum):
    def __init__(self, seqid, start, end):
        Metadatum.__init__(self, 
                           "sequence-region", 
                           "%s %d %d" % (seqid, start, end))
        self.seqid = seqid
        self.start = start
        self.end = end

class FormatError(Exception):
    """Invalid format for GFF file"""
    pass

class Reader:
    """Reads a GFF formatted file"""

    def __init__(self, stream, version="3"):
        self._stream = stream
        self._default_version = version

        # Directives
        self._version = None
        self._references_resolved = True

        # Metadata
        self._metadata = []
        self._comments = []
        self._sequence_regions = []
        self._fasta_string = ""

        # Create record parser table
        self._record_parsers = {"1": self._parse_record_v1,
                                "2": self._parse_record_v2,
                                "2.1": self._parse_record_v2,
                                "2.2": self._parse_record_v2,
                                "2.5": self._parse_record_v2,
                                "3": self._parse_record_v3}

        #  Set default record parser
        if version not in self._record_parsers:
            raise Exception, "Unrecognized GFF default version: " + version
        self._record_parser = self._record_parsers[version]

        # Stage next record so that we read all metadata at top of file
        self._next_rec = None
        self._stage_rec()

    def get_version(self):
        """Returns the format version used for parsing."""
        return self._version or self._default_version

    def is_version_parsed(self):
        """Returns True if the format version was detected in the stream."""
        return self._version is not None

    def get_metadata(self):
        """Returns the metadata read as a list."""
        return self._metadata

    def get_comments(self):
        """Returns comments read as a list."""
        return self._comments

    def get_sequence_regions(self):
        """Returns the sequence region metadata read as a list."""
        return self._sequence_regions

    def get_fasta_string(self):
        """Returns a string of FASTA formatted sequences found at the end of
        the GFF file (version 3 only)"""
        return self._fasta_string

    def are_references_resolved(self):
        """Returns True if record references have all been resolved."""
        return self._references_resolved

    def read(self):
        """Returns the next record or None if there are none left."""
        try:
            return self.next()
        except StopIteration:
            return None

    def read_recs(self, reverse_recs=False):
        """Returns a list of all records that have not yet been read."""
        recs = [rec for rec in self]
        if reverse_recs:
            recs.reverse()
        return recs

    def __iter__(self):
        return self

    def next(self):
        self._stage_rec()
        if self._next_rec is None:
            raise StopIteration
        else:
            rec = self._next_rec
            self._next_rec = None
            return rec

    def _stage_rec(self):
        while self._next_rec is None:
            line = self._stream.readline()

            # Stop when EOF reached
            if line == "":
                return
            # Check for pragma line
            elif line.startswith("##"):
                self._parse_directive(line)
            # Check for comment line
            elif line.startswith("#"):
                self._parse_comment(line)
            # Skip over blank lines
            elif line == "\n":
                pass
            # Check for beginning of FASTA region for v3 formats
            elif line.startswith(">") and self._version == "3":
                self._fasta_string = line + self._stream.read()
            else:
                self._next_rec = self._record_parser(line)

    def _set_record_parser(self):
        try:
            self._record_parser = self._record_parsers[self._version]
        except KeyError:
            self._record_parser = self._record_parsers[self._default_version]
            print >>sys.stderr, "Warning: Unrecognized GFF version (%s). " + \
                "Using default version %s %s." % (self._version, self._default_version)

    def _parse_directive(self, line):
        tokens = line[2:-1].split(None, 1)

        # Skip over empty directive lines
        if not tokens:
            return
            
        # Switch on directive type
        if tokens[0] == "gff-version":
            try:
                self._version = tokens[1]
                self._set_record_parser()
            except IndexError:
                raise FormatError("Invalid gff-version directive: " + line)
        elif tokens[0] == "FASTA":
            # The parser automatically enters FASTA mode with a line starting
            # with '>', so we don't need to do anything here
            pass 
        elif tokens[0] == "#":
            self._references_resolved = True
        # Otherwise this is a metadatum directive
        elif tokens[0] == "sequence-region":
            try:
                seqid, start, end = tokens[1].split()
                seq_region = SequenceRegion(seqid, int(start), int(end))
                self._metadata.append(seq_region)
                self._sequence_regions.append(seq_region)
            except (IndexError, ValueError):
                raise FormatError("Invalid sequence-region directive: " + line)
        else:
            self._metadata.append(Metadatum(*tokens))

    def _parse_comment(self, line):
        # Add line stripped of # prefix and newline to comments list
        self._comments.append(line[1:-1])

    def _parse_record_v1(self, line):
        fields = line[:-1].split('\t', 8)

        if len(fields) == 8:
            attributes = {}
        elif len(fields) == 9:
            attributes = {"group": [fields[8]]}
        else:
            raise FormatError, "Invalid number of fields (should be 8 or 9):\n" + line

        try:
            return GFF(seqid=fields[0],
                          source=fields[1],
                          type=fields[2],
                          start=int(fields[3]),
                          end=int(fields[4]),
                          score=float(fields[5]),
                          strand=parse_maybe_empty(fields[6]),
                          phase=parse_maybe_empty(fields[7], int),
                          attributes=attributes)
        except ValueError, e:
            raise FormatError, "GFF field format error: " + e.message

    def _parse_record_v2(self, line):
        fields = line[:-1].split('\t', 8)

        if len(fields) == 8:
            attributes_string = ""
        elif len(fields) == 9:
            attributes_string = fields[8]
        else:
            raise FormatError, "Invalid number of fields (should be 8 or 9):\n" + line

        try:
            return GFF(seqid=fields[0],
                          source=fields[1],
                          type=fields[2],
                          start=int(fields[3]),
                          end=int(fields[4]),
                          score=parse_maybe_empty(fields[5], float),
                          strand=parse_maybe_empty(fields[6]),
                          phase=parse_maybe_empty(fields[7], int),
                          attributes=self._parse_attributes_v2(attributes_string))
        except ValueError, e:
            raise FormatError, "GFF field format error: " + e.message

    def _parse_record_v3(self, line):
        self._references_resolved = False

        # Strip line
        line = line.strip()
        # Get fields
        fields = line.split('\t')

        if len(fields) != 9:
            raise FormatError, "Invalid number of fields (should be 9):\n" + line

        try:
            return GFF(seqid=url_unquote(fields[0]),
                          source=url_unquote(fields[1]),
                          type=url_unquote(fields[2]),
                          start=int(fields[3]),
                          end=int(fields[4]),
                          score=parse_maybe_empty(fields[5], float),
                          strand=parse_maybe_empty(fields[6]),
                          phase=parse_maybe_empty(fields[7], int),
                          attributes=self._parse_attributes_v3(fields[8]))
        except ValueError, e:
            raise FormatError, "GFF field format error: " + e.message

    def _parse_attributes_v3(self, s):
        attributes = {}

        for pair_string in s.split(";"):
            if (len (pair_string) == 0):
                continue
            try:
                tag, value = pair_string.split("=")
                attributes[url_unquote(tag)] = map(url_unquote,
                                                   value.split(","))
            except ValueError:
                print >>sys.stderr, "WARNING: Invalid attributes string: ", s
#                raise FormatError("Invalid attributes string: " + s)
        return attributes

    def _parse_attributes_v2(self, s):
        attributes = {}
        currentTag = None
        for token in AttributeIterator(s):
            if currentTag is None:
                if isinstance(token, IdentifierToken):
                    currentTag = token.value
                    attributes[currentTag] = []
                else:
                    raise FormatError, "Invalid attributes string: " + s
            elif isinstance(token, SeparatorToken):
                currentTag = None
            elif isinstance(token, CommentToken):
                break
            elif isinstance(token, IdentifierToken):
                attributes[currentTag].append(token.value)
            elif isinstance(token, ValueToken):
                attributes[currentTag].append(token.value)
            else:
                raise FormatError, "Invalid attributes string: " + s
        return attributes

class IdentifierToken:
    pass
class ValueToken:
    pass
class CommentToken:
    pass
class SeparatorToken:
    pass
class UnknownToken:
    pass

class AttributeIterator:
    identifierPat = re.compile(r'\s*([A-Za-z][A-Za-z0-9_]*)')
    freeTextPat = re.compile(r'\s*"(([^"]|(\\"))*)(?<!\\)"')
    valuePat = re.compile(r'\s*([^;# \t\n\r\f\v]+)')
    sepPat = re.compile(r'\s*(;)')
    commentPat = re.compile(r'\s*#(.*)$')

    pats = (identifierPat, freeTextPat, valuePat, sepPat, commentPat)
    tokenClasses = (IdentifierToken, ValueToken, ValueToken,
                    SeparatorToken, CommentToken)

    def __init__(self, s):
        self.s = s.rstrip()
        self.pos = 0
        
    def __iter__(self):
        return self
    
    def next(self):
        if self.pos >= len(self.s):
            raise StopIteration
        
        for (pat, tclass) in zip(AttributeIterator.pats,
                                 AttributeIterator.tokenClasses):
            match = pat.match(self.s, self.pos)
            if match is not None:
                self.pos = match.end(0)
                t = tclass()
                t.value = match.group(1)
                return t
        else:
            return UnknownToken()

def is_integer(x):
    """Returns true if x is of integer type (int or long)."""
    return type(x) in (int, long)

def parse_maybe_empty(s, parse_type=str):
    if s == '.':
        return None
    else:
        return parse_type(s)

def format_maybe_empty(value, empty_str='.'):
    if value is None or value == "":
        return empty_str
    else:
        return str(value)

def quote(s):
    return '"%s"' % str(s)

def url_quote_sub(m):
    return url_quote(m.group(0))

# Added slash '/' to allowable characters
# as well as '(' and ')'
_seqid_pat = re.compile(r'[^a-zA-Z0-9./:^*$@!+_?-|]')
_source_pat = re.compile(r'[^a-zA-Z0-9./\(\): ^*$@!+_?-]')
_type_pat = _source_pat
_tag_pat = re.compile(r'[\t\n\r\f\v;=%&,]')
_value_pat = _tag_pat

class Writer:
    """Writes a GFF formatted file"""

    def __init__(self, stream, version="3", metadata=[]):
        self.stream = stream

        # Create record writer table
        self._record_writers = {"1": self._write_rec_v1,
                                "2": self._write_rec_v2,
                                "2.1": self._write_rec_gtf,
                                "2.2": self._write_rec_gtf,
                                "2.5": self._write_rec_gtf,
                                "3": self._write_rec_v3}

        # Set version
        try:
            self._record_writer = self._record_writers[version]
            self.version = version
            self.write_metadatum(Metadatum("gff-version", version))
        except KeyError:
            raise Exception, "Unrecognized GFF version: " + version

        for metadatum in metadata:
            self.write_metadatum(metadatum)
    
    def write_metadatum(self, metadatum):
        """Writes a metadatum line."""
        if metadatum.value is not None:
            print >>self.stream, "##%s %s" % (metadatum.name, metadatum.value)
        else:
            print >>self.stream, "##%s" % metadatum.name

    def write_comment(self, comment):
        """Writes a comment line."""
        print >>self.stream, "#%s" % comment

    def write(self, rec):
        """Writes a single record."""
        self._record_writer(rec)

    def write_recs(self, recs):
        """Writes a list of records."""
        for rec in recs:
            self.write(rec)

    def _write_rec_v1(self, rec):
        fields = [rec.seqid,
                  rec.source,
                  rec.type,
                  str(rec.start),
                  str(rec.end),
                  format_maybe_empty(rec.score, '0'),
                  format_maybe_empty(rec.strand),
                  format_maybe_empty(rec.phase)]
        if rec.attributes.get('group') is not None:
            fields.append(rec.attributes['group'][0])
        print >>self.stream, '\t'.join(fields)

    def _write_rec_v2(self, rec):
        fields = [rec.seqid,
                  rec.source,
                  rec.type,
                  str(rec.start),
                  str(rec.end),
                  format_maybe_empty(rec.score),
                  format_maybe_empty(rec.strand),
                  format_maybe_empty(rec.phase)]
        if rec.attributes:
            fields.append(self._format_attributes_v2(rec.attributes))
        print >>self.stream, '\t'.join(fields)

    def _write_rec_v3(self, rec):
        _type_pat.sub(url_quote, rec.type)
        fields = [_seqid_pat.sub(url_quote, rec.seqid),
                  _source_pat.sub(url_quote, rec.source),
                  _type_pat.sub(url_quote, rec.type),
                  str(rec.start),
                  str(rec.end),
                  format_maybe_empty(rec.score),
                  format_maybe_empty(rec.strand),
                  format_maybe_empty(rec.phase),
                  format_maybe_empty(self._format_attributes_v3(rec.attributes))]
        print >>self.stream, '\t'.join(fields)

    def _write_rec_gtf(self, rec):
        # The required GTF attributes
        gtf_attributes = ["gene_id", "transcript_id"]
        
        # Make sure that this record has the required GTF attributes
        if not all([attr in rec.attributes for attr in gtf_attributes]):
            gtf_rec = rec.copy()
            for attr in gtf_attributes:
                gtf_rec.attributes.setdefault(attr, [""])
        else:
            gtf_rec = rec

        # GTF is just GFF v2 with some required attributes
        self._write_rec_v2(gtf_rec)
        
    def _format_attributes_v3(self, attributes):
        return ';'.join(["%s=%s" % (_tag_pat.sub(url_quote_sub, tag),
                                    ','.join([_value_pat.sub(url_quote_sub, value)
                                              for value in values]))
                         for tag, values in attributes.items()])

    def _format_attributes_v2(self, attributes):
        return ' '.join([' '.join([tag] + map(quote, values)) + ";"
                          for tag, values in attributes.items()])


def get_inclusive_txn_bounds(gene_hierarchy):
    """
    Get the most inclusive transcription start and end coordinates
    for a given gene.
    """
    mRNA_starts = []
    mRNA_ends = []

    strand = None
    
    for mRNA_id, mRNA_info in gene_hierarchy['mRNAs'].iteritems():
        mRNA_rec = mRNA_info["record"]
        strand = mRNA_rec.strand

        # Get transcription start and end sites
        mRNA_starts.append(mRNA_rec.start)
        mRNA_ends.append(mRNA_rec.end)

    assert(strand != None)

    tx_start = min(mRNA_starts)
    tx_end = max(mRNA_ends)

    # Start must be less than end always
    assert(tx_start < tx_end)
        
    return tx_start, tx_end
