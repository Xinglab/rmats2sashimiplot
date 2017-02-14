import random as pyrand
from scipy import *
from numpy import *
import re
import misopy
from misopy.gff_utils import *
from misopy.parse_csv import *
import pprint


class Interval:
    def __init__(self, start, end):
        self.start = start
        self.end = end
	assert(self.start <= self.end)
	self.len = self.end - self.start + 1
	assert(self.len >= 1)
	
    def __repr__(self):
        return "Interval([%d, %d])" %(self.start, self.end)
    
    def __eq__(self, interval):
	if interval == None: return False
        return self.start == interval.start and self.end == interval.end

    def __ne__(self, interval):
        return not self.__eq(interval)

    def __lt__(self, interval):
        if self.start == interval.start:
            return self.end < interval.end
        return self.start < interval.start

    def contains(self, start, end):
	if self.start <= start and self.end >= end:
	    return True
	return False

    def intersects(self, other):
        if (self.start < other.end
            and self.end > other.start):
            return True
        return False

class Exon(Interval):
    def __init__(self, start, end, label=None, gene=None, seq="",
                 from_gff_record=None):
        Interval.__init__(self, start, end)
	self.gene = gene
	self.label = label
	self.seq = seq
	if self.seq != "":
	    assert(len(self.seq) == len(self.len))
            
        if from_gff_record != None:
            # Load information from a GFF record
            self.load_from_gff_record(from_gff_record)

    def load_from_gff_record(self, gff_record):
        """
        Load exon information from given GFF exon record.
        """
        self.rec = gff_record['record']
        self.parent_rec = gff_record['parent']
        self.start = self.rec.start
        self.end = self.rec.end
        # Use first ID in list
        self.label = self.rec.attributes['ID'][0]

    def __repr__(self):
	gene_label = None
	if self.gene:
	    gene_label = self.gene.label
	return "Exon([%d, %d], id = %s, seq = %s)(ParentGene = %s)" %(self.start, self.end,
								      self.label, self.seq,
                                                                      gene_label)

    def __eq__(self, other):
	if other == None: return False

        # TEST -- removing sequence equality
	if self.start == other.start and self.end == other.end \
	   and self.gene == other.gene:
	    return True
	#if self.seq == other.seq and self.start == other.start and self.end == other.end \
	#   and self.gene == other.gene:
	#    return True
	return False

class Intron(Interval):
    def __init__(self, start, end, label=None, gene=None, seq=""):
        Interval.__init__(self, start, end)
	self.gene = gene
	self.label = label
	self.seq = seq
	if self.seq != "":
	    assert(len(seq) == len(self.len))

    def __repr__(self):
	gene_label = None
	if self.gene:
	    gene_label = self.gene.label
	return "Intron([%d, %d], id = %s)(ParentGene = %s)" %(self.start, self.end,
							      self.label, self.seq,
                                                              self.gene_label)

    def __eq__(self, other):
	if other == None: return False
	if self.seq == other.seq and self.start == other.start and self.end == other.end \
	   and self.gene == other.gene:
	    return True
	return False

class Gene:
    """
    A representation of a gene and its isoforms.  If a gene has two isoforms, make the inclusive
    isoform the first.

    Isoforms are made up of parts, which might be exons or introns.

    isoform_desc is a of lists describing the structure of the isoforms, e.g.

      [['A', 'B', 'A'], ['A', 'A']]

    which creates two isoforms, composed of the 'A' and 'B' parts.
    """    
    def __init__(self, isoform_desc, parts,
                 chrom=None,
                 exons_seq=None,
                 label="",
                 strand="NA",
                 transcript_ids=None):
	self.isoform_desc = isoform_desc
	self.iso_lens = []
	if label == "":
	    self.label = self.get_rand_id(5)
	else:
	    self.label = label
	self.parts = self.create_parts(parts)
	self.isoforms = []
	self.num_parts = len(self.parts)
	self.chrom = chrom
        self.strand = strand
        self.transcript_ids = transcript_ids
        
	# create a set of isoforms
	self.create_isoforms()

        # Use transcript ids if given
        self.assign_transcript_ids()
        
	# The number of exons in each isoform
	self.num_parts_per_isoform = array([iso.num_parts for iso in self.isoforms])

    def isoform_has_part(self, isoform, part):
        """
        Return True if isoform has part, otherwise False.
        """
        for iso_part in isoform.parts:
            if iso_part.start == part.start and \
               iso_part.end == part.end:
                return True
        return False

    def get_const_parts(self):
	"""
	Return all of the gene's constitutive parts, i.e. regions that are shared across all isoforms.
	"""
	const_parts = []
	for part in self.parts:
	    add_part = True
	    for isoform in self.isoforms:
                if not self.isoform_has_part(isoform, part):
                    add_part = False
#		if part not in isoform.parts:
#		    add_part = False
	    if add_part:
		const_parts.append(part)
#	if len(const_parts) == 0:
#	    raise Exception, "Gene %s has no constitutive parts!" %(str(self))
	return const_parts

    def get_alternative_parts(self):
	"""
	Return all of the gene's alternative parts, i.e. non-constitutive regions.
	"""
	const_parts = self.get_const_parts()
	alternative_parts = []
        for part in self.parts:
	    if part not in const_parts:
		alternative_parts.append(part)
	return alternative_parts
            
    def get_rand_id(self, len):
        rand_id = 'G' + "".join([pyrand.choice('abcdefghijklmnopqrstuvwxyz') for n in range(len)])
        return rand_id

    def get_parts_before(self, part):
	"""
	Return all the parts the given part in the gene.
	"""
	parts_before = []
	for p in self.parts:
	    if p == None:
		raise Exception, "Attempting to reference a None part in %s, (part = %s)" \
		      %(str(self), str(part))
	    if p.end < part.start:
		parts_before.append(p)
	    else:
		return parts_before
	return parts_before

    def get_part_number(self, part):
	"""
	Return a part's position (i.e. its number) in the list of parts (0-based).
	"""
	return self.parts.index(part)

    def get_genomic_parts_crossed(self, start, end, read_len=None):
	"""
	Return all parts (by number!) that are crossed in the genomic interval [start, end],
	not including the parts where start and end land.

	If read_len is given, take that into account when considering whether a part was crossed.
	"""
	# find the part where the first coordinate is
	start_part = self.get_part_by_coord(start)
	end_part = self.get_part_by_coord(end)

        ##
        ## NEW: Reads that do not cross any parts might be
        ## intronic reads
        ##
        if start_part == None or end_part == None:
            return []
        
	start_part_num = self.parts.index(start_part)
	end_part_num = self.parts.index(end_part)
	# find parts crossed in between start and end
        parts_crossed = range(start_part_num + 1, end_part_num)	
        if read_len != None:
	    if (end - start) <= read_len:
		return parts_crossed
	return parts_crossed

    def get_part_by_label(self, part_label):
	for part in self.parts:
	    if part_label == part.label:
		return part
	return None

    def part_coords_to_genomic(self, part, part_start, part_end=None):
	"""
	Map the coordinates (start, end) inside part into their corresponding genomic coordinates.
	The end coordinate is optional.
	
	If the given part is the first part in the gene, then these two coordinate systems
	are equivalent.
	"""
	# find parts before the given part and sum their coordinates
	parts_before_len = sum([p.len for p in self.get_parts_before(part)])
	genomic_start = parts_before_len + part_start
	if part_end:
	    genomic_end = parts_before_len + part_end
	    return (genomic_start, genomic_end)
	return genomic_start

    def get_part_by_coord(self, genomic_start):
	"""
	Return the part that contains the given genomic start coordinate.
	"""
	for part in self.parts:
	    if part.start <= genomic_start and genomic_start <= part.end:
		return part
	return None

    def genomic_coords_to_isoform(self, isoform, genomic_start, genomic_end):
	"""
        Get isoform coords and return genomic coords.
	"""
	assert(isoform in self.isoforms)
	# ensure that the parts the coordinates map to are in the isoform
	start_part = self.get_part_by_coord(genomic_start)
	end_part = self.get_part_by_coord(genomic_end)
        
	if start_part == None or end_part == None:
	    return None, None
	# find the isoform coordinates of the corresponding parts
	isoform_start = isoform.part_coord_to_isoform(genomic_start)
	isoform_end = isoform.part_coord_to_isoform(genomic_end)

	return (isoform_start, isoform_end)
		
    def create_parts(self, parts):
	part_counter = 0
	gene_parts = []
	for part in parts:
            gene_part = part

	    gene_part.gene = self
	    gene_parts.append(gene_part)
	    part_counter += 1
	return gene_parts

    def create_isoforms(self):
        self.isoforms = []
	for iso in self.isoform_desc:
            isoform_parts = []
	    isoform_seq = ""
            for part_label in iso:
		# retrieve part 
		part = self.get_part_by_label(part_label)
		if not part:
                    raise Exception, "Invalid description of isoforms: refers to undefined part %s, %s, gene: %s" \
                          %(part, part_label, self.label)
                isoform_parts.append(part)
		isoform_seq += part.seq
	    # make isoform with the given parts
	    isoform = Isoform(self, isoform_parts, seq=isoform_seq)
	    isoform.desc = iso
	    self.isoforms.append(isoform)
	    self.iso_lens.append(isoform.len)
        self.iso_lens = array(self.iso_lens)

    def assign_transcript_ids(self):
        """
        Assign transcript IDs to isoforms.
        """
        if self.transcript_ids != None:
            if len(self.transcript_ids) != len(self.isoforms):
                raise Exception, "Transcript IDs do not match number of isoforms."
            for iso_num, iso in enumerate(self.isoforms):
                curr_iso = self.isoforms[iso_num]
                curr_iso.label = self.transcript_ids[iso_num]
        
    def set_sequence(self, exon_id, seq):
        """
        Set the sequence of the passed in exons to be the given sequence.
        """
	return Exception, "Unimplemented method."

    def align_read_pair_with_cigar(self, left_cigar, genomic_left_read_start,
                                   genomic_left_read_end, right_cigar,
                                   genomic_right_read_start,
                                   genomic_right_read_end, read_len, 
                                   overhang=1):

        alignment = []
        iso_frag_lens = []
        
        left = self.align_read_to_isoforms_with_cigar(
            left_cigar, genomic_left_read_start, genomic_left_read_end,
            read_len, overhang)
        right = self.align_read_to_isoforms_with_cigar(
            right_cigar, genomic_right_read_start, genomic_right_read_end,
            read_len, overhang)
    
        for lal, lco, ral, rco in zip(left[0], left[1], right[0], right[1]):
            if lal and ral:
                alignment.append(1)
                iso_frag_lens.append(rco[1]-lco[0]+1)
            else:
                alignment.append(0)
                iso_frag_lens.append(-Inf)

        return (alignment, iso_frag_lens)

#     def align_read_pair(self, genomic_left_read_start, genomic_left_read_end,
# 			genomic_right_read_start, genomic_right_read_end, overhang=1,
# 			read_len=36):
# 	"""
# 	Align a paired-end read to all of the gene's isoforms.
	
# 	Return an alignment binary vector of length K, where K is the number of isoforms, as well as
# 	a vector with the fragment lengths that correspond to each alignment to an isoform.

# 	When a read does not align to a particular isoform, the fragment length for that alignment is
# 	denoted with -Inf.
# 	"""
# 	# align each of the pairs independently to all isoforms, and then compute the fragment
# 	# lengths that correspond to each alignment
# 	alignment = []
# 	iso_frag_lens = []

# #        print "Aligning: ", (genomic_left_read_start, genomic_left_read_end), \
# #              " - ", (genomic_right_read_start, genomic_right_read_end)
        
# 	# align left read
# 	(left_alignment, left_isoform_coords) = self.align_read_to_isoforms(genomic_left_read_start, genomic_left_read_end,
# 									    overhang=overhang, read_len=read_len)
# 	# align right read
# 	(right_alignment, right_isoform_coords) = self.align_read_to_isoforms(genomic_right_read_start,
# 									      genomic_right_read_end,
# 									      overhang=overhang, read_len=read_len)

# 	num_isoforms = len(self.isoforms)
# 	for n in range(num_isoforms):
# 	    # check that both reads align to the isoform with no overhang violation
# 	    if not (left_alignment[n] == right_alignment[n] and right_alignment[n] != 0):
# 		alignment.append(0)
# 		iso_frag_lens.append(-Inf)
# 		continue
# #            else:
# #                print "Not both reads align to the same isoform"
# 	    # compute fragment length for each isoform, which is the isoform start coordinate of
# 	    # the right read minus the isoform start coordinate of the left read
# 	    frag_len = right_isoform_coords[n][1] - left_isoform_coords[n][0] + 1
# 	    if frag_len < 0:
# 		# negative fragment length		
# 		print "Warning: Negative fragment length during alignment of ", genomic_left_read_start, \
# 		      genomic_right_read_start, " to isoform: ", self.isoforms[n]
# 		alignment.append(0)
# 		iso_frag_lens.append(-Inf)
# 		continue
# 	    # both reads align with no overhang violations and the fragment length is reasonable
# 	    alignment.append(1)
# 	    iso_frag_lens.append(frag_len)
# 	return (alignment, iso_frag_lens)

    def align_reads_to_isoforms(self, read_genomic_coords, overhang=1, read_len=36):
	alignments = []
	isoforms_coords = []
	for read_coords in read_genomic_coords:
	    genomic_read_start, genomic_read_end = read_coords
	    alignment, isoform_coords = self.align_read_to_isoforms(genomic_read_start, genomic_read_end,
                                                                    overhang=overhang, read_len=read_len)
	    alignments.append(alignment)
	    isoforms_coords.append(isoform_coords)
	return (array(alignments), isoforms_coords)

    def align_read_to_isoforms_with_cigar(self, cigar, genomic_read_start,
                                          genomic_read_end, read_len, overhang_len):
        """
        Align a single-end read to all of the gene's isoforms.
        Use the cigar string of the read to determine whether an
        isoform matches
        """
        alignment = []
        isoform_coords = []
        for isoform in self.isoforms:
            iso_read_start, iso_read_end = self.genomic_coords_to_isoform(isoform,
                                                                          genomic_read_start,
                                                                          genomic_read_end)
            isocigar = isoform.get_local_cigar(genomic_read_start, read_len)

            # Check that read is consistent with isoform and that the overhang
            # constraint is met
            if (isocigar and isocigar == cigar) and \
                   isoform.cigar_overhang_met(isocigar, overhang_len):
                alignment.append(1)
                isoform_coords.append((iso_read_start, iso_read_end))
            else:
                alignment.append(0)
                isoform_coords.append(None)
                
        return (alignment, isoform_coords)        


#     def align_read_to_isoforms(self, genomic_read_start, genomic_read_end, overhang=1, read_len=36):
# 	"""
# 	Align a single-end read to all of the gene's isoforms.

# 	Return an alignment as well as a set of isoform coordinates, for each isoform, corresponding
# 	to the places where the read aligned.

#         When the read violates overhang or the read doesn't align to a particular
# 	isoform, the coordinate is set to None.
# 	"""
# 	alignment = []
# 	isoform_coords = []
# 	genomic_parts_crossed = self.get_genomic_parts_crossed(genomic_read_start,
#                                                                genomic_read_end,
# 							       read_len=read_len)

#         if len(genomic_parts_crossed) == 0:
#             print "zero genomic parts crossed"

#         ##
#         ## NEW: If no genomic parts are crossed, must be intronic read
#         ##
#         # if len(genomic_parts_crossed) == 0:
#         #     alignment = [0] * len(self.isoforms)
#         #     isoform_coords = [None] * len(self.isoforms)
#         #     return (alignment, isoform_coords)
        
# 	for isoform in self.isoforms:
# 	    # check that parts aligned to in genomic coordinates exist
# 	    # in the current isoform
# 	    iso_read_start, iso_read_end = self.genomic_coords_to_isoform(isoform,
# 									  genomic_read_start,
# 									  genomic_read_end)
# 	    if iso_read_start == None or iso_read_end == None:
# 		alignment.append(0)
# 		isoform_coords.append(None)
# 		continue
# 	    # genomic parts have matching parts in isoform. Now check that
# 	    # that they cross the same junctions
# 	    iso_parts_crossed = isoform.get_isoform_parts_crossed(iso_read_start, iso_read_end)

# 	    if iso_parts_crossed != genomic_parts_crossed:
# 		alignment.append(0)
# 		isoform_coords.append(None)
# 		continue
            
# 	    # check that overhang violation is met on outer parts crossed (as long as
# 	    # parts that are crossed in between are greater or equal to overhang constraint,
# 	    # no need to check those)
# 	    if overhang > 1:
# 		start_part, start_part_coord = isoform.get_part_by_coord(iso_read_start)
#                 if (start_part.end - (start_part_coord + start_part.start) + 1) < overhang:
# 		    alignment.append(0)
# 		    isoform_coords.append(None)		 
# 		    continue
# 		end_part, end_part_coord = isoform.get_part_by_coord(iso_read_end)
# 		if ((end_part_coord + end_part.start) - end_part.start) + 1 < overhang:
# 		    alignment.append(0)
# 		    isoform_coords.append(None)		    
# 		    continue
# 	    # overhang is met and read aligns to the isoform
# 	    alignment.append(1)
# 	    # register coordinates
# 	    isoform_coords.append((iso_read_start, iso_read_end))
# 	return (alignment, isoform_coords)

    def align_read(self, genomic_read_start, genomic_read_end, overhang=1, read_len=36):
	"""
	Align a single-end read to all of the gene's isoforms.
	Return an alignment binary vector of length K, where K is the number of isoforms.
	The ith position in this vector has a 1 if the read aligns to the given isoform,
	and 0 otherwise.

	A read is of the form:

	  (genomic_start_coord, genomic_end_coord)
	"""
	alignment = []
	# align all the reads to isoforms and return an alignment back, as well as the set of
	# genomic coordinates for each isoform that the read aligns to.
	alignment, aligned_genomic_coords = self.align_read_to_isoforms(genomic_read_start, genomic_read_end,
									overhang=overhang)
	# get the parts that are crossed in genomic coordinates space, taking into account
	# the read's length
	category = None
	two_iso_alignment = None
	# Deal with the two-isoform special case: categorize reads into NI, NE, NB
	if len(self.isoforms) == 2:
	    if alignment == [1, 0]:
		# NI
		two_iso_alignment = [1, 0, 0]
		# Check if it aligns to upstream or downstream inclusion junction
		# We know that overhang violation hasn't been made here based on the alignment
		# first convert the genomic coordinates of read to the first isoform's coordinates
		isoform1 = self.isoforms[0]
		# find isoform coordinate of genomic read start
		iso1_read_start, c1 = self.genomic_coords_to_isoform(isoform1, genomic_read_start,
								     genomic_read_start)
		# find isoform coordinate of genomic read end
		iso1_read_end, c2 = self.genomic_coords_to_isoform(isoform1, genomic_read_end,
								   genomic_read_end)
		# find which parts these coordinates land in
		iso1_read_start_part, c1 = isoform1.get_part_by_coord(iso1_read_start)
		iso1_read_end_part, c2 = isoform1.get_part_by_coord(iso1_read_end)
		# if the read starts in the first part of the isoform and
		# ends in the second part of the isoform, then it's an upstream
		# inclusion junction read
		if iso1_read_start_part == isoform1.parts[0] and iso1_read_end_part == isoform1.parts[1]:
		    category = 'upincjxn'
		elif iso1_read_start_part == isoform1.parts[1] and iso1_read_end_part == isoform1.parts[2]:
		    # if the read starts in the second part of the isoform and
		    # ends in the third, then it's a downstream inclusion
		    # junction read
		    category = 'dnincjxn'
		elif iso1_read_start_part == isoform1.parts[1] and iso1_read_end_part == isoform1.parts[1]:
		    # if the read starts and ends in the skipped exon body,
		    # then it's a body read
		    category = 'body'
		else:
		    # If the read is not in one of those categories, the isoform must have more than three parts
		    assert(len(isoform1.parts) > 3)
		# if the read doesn't fall into either of these categories, it can't possibly be
		# an inclusion read
#		if category == None:
#		    raise Exception, "Incoherent inclusion read: not upincjxn, dnincjxn, or body! %s" \
#			  %(str(iso1_read_start) + ' - ' + str(iso1_read_end))
	    elif alignment == [0, 1]:
		# NE
		two_iso_alignment = [0, 1, 0]
	    elif alignment == [1, 1]:
		# NB
		two_iso_alignment = [0, 0, 1]
	    else:
		# Overhang violation
		two_iso_alignment = [0, 0, 0]
	    return (two_iso_alignment, category)
	return (alignment, category)

    def align_paired_end_reads(self, reads, overhang=1):
	"""
	Take a set of paired-end reads parameterized by their genomic coordinates
	and align them to all of the gene's isoforms.

	For each read, return a pair where the first element is the alignment to all the
	isoforms (a binary vector, with 1 in the ith position of the read aligns to the
	ith isoform and 0 otherwise) and the second element is the length of the
	fragment lengths that correspond to each alignment (-Inf if the read does not
	align to the given isoform.)
	"""
	aligned_reads = []
	for read in reads:
	    # align read to all of the isoforms
	    (alignment, isoform_coords) = self.align_read_pair(read[0], read[1], read[2], read[3], overhang=overhang)
	    frag_lens = [c2 - c1 + 1 for c1, c2 in isoform_coords]
	    aligned_reads.append(array([alignment, frag_lens]))
	return aligned_reads
	
    def align(self, seq, overhang=None):
        """
        Given a short sequence, return a list of size len(self.isoforms) that says
        for each isoform if the sequence is a substring of it (denoted 1) or not
        (denoted 0).
        """
        # if no overhang constraints are given, align without regard for overhang violations
        alignment = []
        if not overhang:
            for iso in self.isoforms:
                alignment.append(1 if seq in iso.seq else 0)
            return alignment
        category = None
        # take overhang into account
        for iso in self.isoforms:
            # find read starting position in the isoform
            read_start = iso.seq.find(seq)
	    if read_start == -1:
                alignment.append(0)
                continue
            split_iso = iso.seq[read_start:]
            prev_part = 0
            oh_viol = False
	    for part in iso.parts:
                # find beginning of next part
                next_part = split_iso.find(part.seq)
                if next_part == -1:
                    continue
                # check that there is no overhang violation
                part_segment = seq[prev_part:next_part]
                remain_seq = seq[next_part:]
                prev_part = next_part
                if len(part_segment) != 0 and len(part_segment) < 4:
                    oh_viol = True
                    alignment.append(0)
                    break
            if not oh_viol and remain_seq != '':
                if len(remain_seq) < 4:
                    oh_viol = True
                    alignment.append(0)
            if not oh_viol:
                alignment.append(1)
        # Deal with the two isoform case. Classify each read into a category 
        # and then return the counts [ni, ne, nb].
        if len(self.isoforms) == 2:
            if re.match("^0000.*1111.*$", seq) != None:
                category = 'upincjxn'
            elif re.match("^11111*$", seq) != None:
                category = 'body'
            elif re.match("^11111*22222*$", seq) != None:
                category = 'dnincjxn'
            if alignment == [1, 0]:
                return ([1, 0, 0], category)
            elif alignment == [0, 1]:
                return ([0, 1, 0], category)
            elif alignment == [1, 1]:
                return ([0, 0, 1], category)
            elif alignment == [0, 0]:
                # overhang violation
                return ([0, 0, 0], category)                
        return (alignment, category)
        
    def get_iso(self, isoform_num):
        return self.isoforms[isoform_num]

    def avg_iso_len(self):
        """
        Return the gene's average isoform length.
        """
        iso_lens = [i['len'] for i in self.isoforms]
        return mean(iso_lens)

    def __str__(self):
        return "gene_id: %s\nisoforms: %s" %(self.label, self.isoforms)

    def __repr__(self):
        return self.__str__()  

class Isoform:
    def __init__(self, gene, parts,
                 seq=None,
                 label=None):
	"""
	Builds an isoform given an isoform description.
	"""
	self.gene = gene
	# ordered list of isoform parts (exons and introns)
	self.parts = parts
	self.num_parts = len(parts)
	self.len = sum([part.len for part in parts])
	self.seq = seq
        self.label = label
	# the genomic coordinates of the isoform are defined as the start coordinate
	# of the first part and the last coordinate of the last part
	first_part = self.parts[0]
	last_part = self.parts[-1]
	self.genomic_start = first_part.start
	self.genomic_end = last_part.end
	
    def get_parts_before(self, part):
	"""
	Return all the parts the given part in the isoform:
	"""
	parts_before = []
	for p in self.parts:
	    if p.end < part.start:
		parts_before.append(p)
	    else:
		return parts_before
	return parts_before

    def get_part_by_coord(self, start_coord):
	"""
	Get the part that the *given isoform start coordinate* lands in, and the corresponding
	part-based coordinate.
	"""
	isoform_interval_start = 0
	isoform_interval_end = 0
	prev_part = None
	for part in self.parts:
	    # count the parts in isoform coordinate space and see in which part
	    # the given coordinate falls
	    isoform_interval_end += part.len - 1
	    # check that the part contains the single point
	    if isoform_interval_start <= start_coord and start_coord <= isoform_interval_end:
		# find parts before and sum them up to get the part-based coordinate
		# the part based coordinate is start_coord - previous part lengths
		prev_parts_sum = sum([p.len for p in self.get_parts_before(part)])
		part_start = start_coord - prev_parts_sum
		return part, part_start
	    # add one to move to next part
	    isoform_interval_end += 1
	    isoform_interval_start = isoform_interval_end
	return (None, None)

    def get_isoform_parts_crossed(self, start, end):
	"""
	Return all parts (by number!) that are crossed in the isoform interval [start, end],
	not including the parts where start and end land.
	"""
	# find the part where the first coordinate is
	start_part, s1 = self.get_part_by_coord(start)
	end_part, s2 = self.get_part_by_coord(end)
	start_part_num = self.parts.index(start_part)
	end_part_num = self.parts.index(end_part)
	# find parts crossed in between start and end	
	return range(start_part_num + 1, end_part_num)

    def cigar_overhang_met(self, cigar, overhang_len):
        """
        Check that the overhang constraint is met in each
        match condition of the read.
        """
        overhang_met = True
        for c in cigar:
            # If it's the match (M) part of the cigar
            # and the match length is less than the overhang
            # constraint, then the constraint is violated
            if (c[0] == 0) and (c[1] < overhang_len):
                return False
        return overhang_met

    def get_local_cigar(self, start, read_len):
        """
        Calculate a CIGAR string for a hypothetical read at a given start position, with a given read length"""
        # If the read starts before or after the isoform, then it does not fit
        if start < self.parts[0].start or self.parts[-1].end < start:
            return None
        # Look for the exon where the read starts
        found = None
        for i, p in enumerate(self.parts):
            if p.start <= start and start <= p.end:
                found = i
                break
        if found == None:
            return None
        
        # Create CIGAR string
        cigar = []
        rl = read_len
        st = start
        for i in range(found, len(self.parts)):
            # the rest is on this exon?
            if rl <= self.parts[i].end - st + 1:
                cigar.append((0, rl))
                return cigar
            # the next exon is needed as well
            else:
                # is there a next exon?
                if i+1 == len(self.parts):
                    return None
                cigar.append((0, self.parts[i].end - st + 1))
                cigar.append((3, self.parts[i+1].start - self.parts[i].end - 1))
                rl = rl - (self.parts[i].end - st + 1)
                st = self.parts[i+1].start
        return cigar

    def part_coord_to_isoform(self, part_start):
	"""
	Get the isoform coordinate that the *given part_start coordinate* lands in.
	"""
	isoform_interval_start = 0
	isoform_coord = None
	for part in self.parts:
            if part.contains(part_start, part_start):
		isoform_coord = isoform_interval_start + (part_start - part.start)
		return isoform_coord
	    isoform_interval_start += part.len
	return isoform_coord

    def isoform_coords_to_genomic(self, isoform_start, isoform_end):
	"""
	Map coordinates of isoform to genomic coordinates.
	"""
	# get the part that each coordinate point lands in
	start_part, start_part_coord = self.get_part_by_coord(isoform_start)
	end_part, end_part_coord = self.get_part_by_coord(isoform_end)
        
	# retrieve the corresponding genomic coordinates
	genomic_start = self.gene.part_coords_to_genomic(start_part, start_part_coord)
	genomic_end = self.gene.part_coords_to_genomic(end_part, end_part_coord)
	return (genomic_start, genomic_end)

    def __repr__(self):
	parts_str = str([p.label for p in self.parts])
	return "Isoform(gene = %s, g_start = %d, g_end = %d, len = %d,\n parts = %s)" \
	       %(self.gene.label, self.genomic_start, self.genomic_end, self.len, parts_str)


def pretty(d, indent=0):
   for key, value in d.iteritems():
      print '  ' * indent + str(key)
      if isinstance(value, dict):
         pretty(value, indent+1)
      else:
         print '  ' * (indent+1) + str(value)
        

def printTree(tree, depth = 0):
    if tree == None or not type(tree) == dict:
        print "\t" * depth, tree
    else:
        for key, val in tree.items():
            print "\t" * depth, key
            printTree(val, depth+1)
            

def print_gene_hierarchy(gene_hierarchy):
    pretty(gene_hierarchy)
#    pp = pprint.PrettyPrinter(indent=4)
#    pp.pprint(gene_hierarchy)

def load_genes_from_gff(gff_filename,
                        include_introns=False,
                        reverse_recs=False,
                        suppress_warnings=False):
    """
    Load all records for a set of genes from a given GFF file.
    Parse each gene into a Gene object.
    """
    gff_db = GFFDatabase(gff_filename,
                         include_introns=include_introns,
                         reverse_recs=reverse_recs)
    # dictionary mapping gene IDs to the list of all their relevant records
    gff_genes = {}

    num_genes = 0

    for gene in gff_db.genes:
	gene_records, gene_hierarchy = gff_db.get_genes_records([gene.get_id()])

        # Record the gene's GFF record
        gene_label = gene.get_id()

        if gene_label not in gene_hierarchy:
            if not suppress_warnings:
                print "Skipping gene %s..." %(gene_label)
            continue
        
        gene_hierarchy[gene_label]['gene'] = gene

        # Make a gene object out of the GFF records
        gene_obj = make_gene_from_gff_records(gene_label,
                                              gene_hierarchy[gene_label],
                                              gene_records)
        if gene_obj == None:
            if not suppress_warnings:
                print "Cannot make gene out of %s" %(gene_label)
            continue
        gff_genes[gene.get_id()] = {'gene_object': gene_obj,
                                    'hierarchy': gene_hierarchy}

        if (num_genes % 5000) == 0:
            if not suppress_warnings:
                print "Through %d genes..." %(num_genes)
        num_genes += 1

    num_genes = len(gff_genes)
    if not suppress_warnings:
        print "Loaded %d genes" %(num_genes)

    return gff_genes


def make_gene_from_gff_records(gene_label,
                               gene_hierarchy,
                               gene_records):
    """
    Make a gene from a gene hierarchy.
    """
    mRNAs = gene_hierarchy['mRNAs']

    # Each transcript is a set of exons
    transcripts = []
    isoform_desc = []

    chrom = None
    strand = "NA"

    # Iterate through mRNAs in the order in which they were given in the
    # GFF file
    transcript_ids = [rec.get_id() for rec in gene_records \
                      if (rec.type == "mRNA" or rec.type == "transcript")]

    if len(transcript_ids) == 0:
        raise Exception, "Error: %s has no transcripts..." \
              %(gene_label)

    num_transcripts_with_exons = 0

    used_transcript_ids = []

    for transcript_id in transcript_ids:
        transcript_info = mRNAs[transcript_id]
        transcript_rec = transcript_info['record']

        chrom = transcript_rec.seqid
        strand = transcript_rec.strand
        transcript_exons = transcript_info['exons']
        exons = []

        if len(transcript_exons) == 0:
            print "%s has no exons" %(transcript_id)
            continue

        # Record how many transcripts we have with exons children
        # (i.e., usable transcripts)
        num_transcripts_with_exons += 1

        for exon_id, exon_info in transcript_exons.iteritems():
            exon_rec = exon_info['record']

            exon = Exon(exon_rec.start, exon_rec.end, from_gff_record={'record':
                                                                       exon_rec,
                                                                       'parent':
                                                                       transcript_rec})
            exons.append(exon)

        # Sort exons by their start coordinate
        exons = sorted(exons, key=lambda e: e.start)

        # Exons that make up a transcript
        transcripts.append(exons)

        # Get exon labels to make transcript's description
        exon_labels = [exon.label for exon in exons]

        # Delimiter for internal representation of isoforms
        #iso_delim = "_"
        
        # The transcript's description
        #for label in exon_labels:
        #    if iso_delim in label:
        #        raise Exception, "Cannot use %s in naming exons (%s) in GFF." %(iso_delim,
        #                                                                        label)
            
        #desc = iso_delim.join(exon_labels)
        isoform_desc.append(exon_labels)
        # Record transcript ids that are not skipped
        used_transcript_ids.append(transcript_id)

    #if num_transcripts_with_exons < 2:
    #    print "WARNING: %s does not have at least two mRNA/transcript entries " \
    #          "with exons. Skipping over..." %(gene_label)
    #    return None

    # Compile all exons used in all transcripts
    all_exons = []
    [all_exons.extend(transcript) for transcript in transcripts]

    # Prefix chromosome with "chr" if it does not have it already
    #if not chrom.startswith("chr"):
    #    chrom = "chr%s" %(chrom)

    gene = Gene(isoform_desc, all_exons,
                label=gene_label,
                chrom=chrom,
                strand=strand,
                transcript_ids=used_transcript_ids)

    return gene
                        

def make_gene(parts_lens, isoforms, chrom=None):
    """
    Make a gene out of the given parts lengths, where isoforms are a
    list of list of numbers, where each list of numbers is an isoform
    (the numbers corresponding to the exon parts in part_lens --
    *one-based* index).
    """
    parts = []
    start_genomic = 0
    part_num = 1
    for part_len in parts_lens:
	end_genomic = start_genomic + part_len - 1
	exon = Exon(start_genomic, end_genomic, label=str(part_num))
	parts.append(exon)
	part_num += 1
	start_genomic = end_genomic + 1
    isoform_desc = []
    for iso in isoforms:
	desc = "_".join([str(iso_name) for iso_name in iso])
	isoform_desc.append(desc)
    gene = Gene(isoform_desc, parts, chrom=chrom)
    return gene

def se_event_to_gene(up_len, se_len, dn_len, chrom,
                     label=None):
    """
    Parse an SE event to a gene structure.
    """
    exon1_start = 0
    exon1_end = up_len - 1
    
    exon2_start = exon1_end + 1
    exon2_end = exon2_start + (se_len - 1)
    
    exon3_start = exon2_end + 1
    exon3_end = exon3_start + (dn_len - 1)
    up_exon = Exon(exon1_start, exon1_end, label='A')
    se_exon = Exon(exon2_start, exon2_end, label='B')
    dn_exon = Exon(exon3_start, exon3_end, label='C')
    parts = [up_exon, se_exon, dn_exon]
    gene = Gene([['A', 'B', 'C'], ['A', 'C']], parts, label=label,
                chrom=chrom)
    return gene

def tandem_utr_event_to_gene(core_len, ext_len, chrom, label=None):
    """
    Parse a tandem UTR event to a gene structure.
    """
    exon1_start = 0
    exon1_end = core_len - 1

    exon2_start = exon1_end + 1
    exon2_end = exon2_start + (ext_len - 1)

    core_exon = Exon(exon1_start, exon1_end, label='TandemUTRCore')
    ext_exon = Exon(exon2_start, exon2_end, label='TandemUTRExt')
    parts = [core_exon, ext_exon]
    gene = Gene([['TandemUTRCore', 'TandemUTRExt'], ['TandemUTRCore']], parts,
                label=label, chrom=chrom)
    return gene

def make_proximal_distal_exon_pair(proximal_exons, distal_exons):
    """
    Make one large distal exon out of a small 
    """
    return

def afe_ale_event_to_gene(proximal_exons, distal_exons, event_type,
                          chrom, read_len=None, overhang_len=None,
                          label=None):
    """
    Parse an AFE/ALE event to a gene.
    """
    # Extend each exon by the junction that can be made between it
    # and the gene body, based on read length and overhang
    if read_len != None and overhang_len != None:
        #num_junction_positions = read_len - (2 * (overhang_len - 1))
        num_junction_positions = read_len
    else:
        num_junction_positions = 0

    # Distal exon - farther away from the gene body
    distal_exon_start = 0
    sum_distal_exons_lens = sum([distal_exon['len'] for distal_exon \
                                 in distal_exons])
    sum_distal_exons_lens += num_junction_positions
    distal_exon_end = sum_distal_exons_lens - 1
    distal_exon = Exon(distal_exon_start, distal_exon_end,
                       label='%sDistal' %(event_type))

    # Proximal exon - closer to the gene body
    proximal_exon_start = distal_exon_end + 1
    sum_proximal_exons_lens = sum([proximal_exon['len'] for proximal_exon \
                                   in proximal_exons])
    sum_proximal_exons_lens += num_junction_positions
    proximal_exon_end = proximal_exon_start + (sum_proximal_exons_lens - 1)
    proximal_exon = Exon(proximal_exon_start, proximal_exon_end,
                         label='%sProximal' %(event_type))
    parts = None

    if event_type == 'AFE':
        parts = [distal_exon, proximal_exon]
    else:
        raise Exception, "Parsing wrong event type, %s" %(event_type)

    # Make it so proximal isoform is always first
    gene = Gene(['%sProximal' %(event_type), '%sDistal' %(event_type)],
                parts, chrom=chrom, label=label)
    return gene

if __name__ == '__main__':
    pass
    
