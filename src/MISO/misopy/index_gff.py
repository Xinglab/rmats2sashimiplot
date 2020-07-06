# -*- mode: python; -*-
##
## Script to build an indexed representation of a GFF file for efficient
## retrieval of genes
##
import os
import sys
import time
import glob
import shelve

from collections import defaultdict
# Add misopy path
miso_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, miso_path)


import misopy
import misopy.gff_utils as gff_utils
import misopy.pickle_utils as pickle_utils
import misopy.Gene as gene_utils
import misopy.misc_utils as misc_utils
print(misopy)


COMPRESS_PREFIX = misc_utils.COMPRESS_PREFIX

def compress_event_name(event_name,
                        prefix=COMPRESS_PREFIX):
    event_hash = hash(event_name)
    compressed_event_name = "%s_%s" %(prefix, event_hash)
    return compressed_event_name


def serialize_genes(gff_genes,
                    gff_filename,
                    output_dir,
                    compress_id=False):
    """
    Output genes into pickle files by chromosome, by gene.

    If asked, use compressed IDs (hashes) of the 'ID=' field in the GFF.
    """
    genes_by_chrom = defaultdict(dict)

    # Split up genes by chromosome
    for gene_id, gene_info in gff_genes.iteritems():
        gene_obj = gene_info["gene_object"]
        gene_hierarchy = gene_info["hierarchy"]
        genes_by_chrom[gene_obj.chrom][gene_id] = \
            {'gene_object': gene_obj,
             'hierarchy': gene_hierarchy}
        if compress_id:
            gene_compressed_id = compress_event_name(gene_id)
            # Store compressed ID
            genes_by_chrom[gene_obj.chrom][gene_id]['compressed_id'] \
                = gene_compressed_id

    # Mapping from gene IDs to pickled filename
    gene_id_to_filename = {}
    # Mapping from compressed IDs (hashes) to gene IDs
    compressed_id_to_gene_id = {}

    # Serialize all the genes in each chromosome into their
    # own directory
    for chrom, chrom_genes in genes_by_chrom.iteritems():
        if chrom.startswith("chr"):
            chrom_dir_name = chrom
        else:
            # Add chr-prefix for ease of finding directory
            # in downstream steps.
            chrom_dir_name = "chr%s" %(str(chrom))

        # Make directory for chromosome if it doesn't already exist
        chrom_dir = os.path.join(output_dir, chrom_dir_name)
        if not os.path.isdir(chrom_dir):
            print "Making directory: %s" %(chrom_dir)
            os.makedirs(chrom_dir)

        t1 = time.time()
        # Serialize each gene into a separate file
        num_genes = len(genes_by_chrom[chrom])

        for gene_id, gene_info in genes_by_chrom[chrom].iteritems():
            gene_compressed_id = None
            if compress_id:
                gene_compressed_id = \
                    genes_by_chrom[chrom][gene_id]['compressed_id']
                gene_filename = \
                    os.path.abspath(os.path.join(chrom_dir,
                                                 "%s.pickle" \
                                                 %(gene_compressed_id)))
            else:
                gene_filename = \
                    os.path.abspath(os.path.join(chrom_dir,
                                                 "%s.pickle" %(gene_id)))
            # Write each gene/event's pickle file
            pickle_utils.write_pickled_file({gene_id:
                                             genes_by_chrom[chrom][gene_id]},
                                            gene_filename)
            # Record what filename was associated with this gene ID
            gene_id_to_filename[gene_id] = gene_filename
            # Record compressed ID (hash) to gene ID
            if gene_compressed_id is not None:
                compressed_id_to_gene_id[gene_compressed_id] = gene_id

        t2 = time.time()
        print "  - Chromosome serialization took %.2f seconds" %(t2 - t1)

    # Shelve the mapping from gene ids to filenames
    shelved_filename = os.path.join(output_dir,
                                    "genes_to_filenames.shelve")
    shelved_data = shelve.open(shelved_filename)
    for k, v in gene_id_to_filename.iteritems():
        shelved_data[k] = v
    shelved_data.close()

    # Shelve the mapping from compressed gene ids to gene ids
    shelved_filename = os.path.join(output_dir,
                                    "compressed_ids_to_genes.shelve")
    shelved_data = shelve.open(shelved_filename)
    for k, v in compressed_id_to_gene_id.iteritems():
        shelved_data[k] = v
    shelved_data.close()

    # Output a list of genes in ordinary GFF format
    genes_filename = os.path.join(output_dir, "genes.gff")
    print "Outputting gene records in GFF format..."
    print "  - Output file: %s" %(genes_filename)
    with open(gff_filename) as gff_in:
        with open(genes_filename, "w") as gff_out:
            for line in gff_in:
                if line.startswith("#"): continue
                record_type = line.strip().split("\t")[2]
                if record_type == "gene":
                    gff_out.write(line)


def index_gff(gff_filename, output_dir,
              compress_id=False):
    """
    Index the given GFF and placed the indexed representation
    in the output directory.
    """
    print "Indexing GFF..."
    if compress_id:
        print "  - Using compressed IDs to create indexed filenames."
    # First check that the GFF is not already indexed
    indexed_files = glob.glob(os.path.join(output_dir, "chr*"))
    if len(indexed_files) >= 1:
        print "%s appears to already be indexed. Aborting." \
            %(gff_filename)
        return

    print "  - GFF: %s" %(gff_filename)
    print "  - Outputting to: %s" %(output_dir)
    overall_t1 = time.time()
    t1 = time.time()
    gff_genes = gene_utils.load_genes_from_gff(gff_filename)
    t2 = time.time()
    print "  - Loading of genes from GFF took %.2f seconds" %(t2 - t1)

    t1 = time.time()
    serialize_genes(gff_genes,
                    gff_filename,
                    output_dir,
                    compress_id=compress_id)
    t2 = time.time()
    print "  - Serialization of genes from GFF took %.2f seconds" %(t2 - t1)
    overall_t2 = time.time()
    print "Indexing of GFF took %.2f seconds." %(overall_t2 - overall_t1)


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--index", dest="index_gff", nargs=2, default=None,
                      help="Index the given GFF. Takes as arguments as GFF filename "
                      "and an output directory.")
    parser.add_option("--compress-id", dest="compress_id", default=False,
                      action="store_true",
                      help="Use the compressed version of the GFF \'ID=\' "
                      "field rather than the ID itself when creating "
                      ".miso output filenames.")
    (options, args) = parser.parse_args()

    if options.index_gff != None:
        gff_filename = \
            os.path.abspath(os.path.expanduser(options.index_gff[0]))
        output_dir = \
            os.path.abspath(os.path.expanduser(options.index_gff[1]))
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        index_gff(gff_filename, output_dir,
                  compress_id=options.compress_id)
    else:
        print "Indexer of GFF files for use with MISO."
        print "Need to pass --index, for example:\n"
        print "index_gff --index annotation.gff indexed_annotation/"


if __name__ == '__main__':
    main()
