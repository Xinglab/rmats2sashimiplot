import os

import misopy
import misopy.gff_utils as gff_utils

def parseGene(pickle_filename, event):
    """
    Parse a pickled gene.
    """
    if not os.path.isfile(pickle_filename):
        raise Exception, "Error: no filename %s" %(pickle_filename)
    gff_genes = gff_utils.load_indexed_gff_file(pickle_filename)

    if gff_genes == None:
        raise Exception, "Error: could not load genes from %s" \
              %(pickle_filename)

    exon_starts = []
    exon_ends = []
    mRNAs = []
    chrom = None
    for gene_id, gene_info in gff_genes.iteritems():
        if event == gene_id:
            gene_obj = gene_info['gene_object']
            gene_hierarchy = gene_info['hierarchy']
            tx_start, tx_end = gff_utils.get_inclusive_txn_bounds(\
                gene_hierarchy[gene_id])
            chrom = gene_obj.chrom

            for mRNA_id, mRNA_info in gene_hierarchy[gene_id]['mRNAs'].iteritems():
                mRNA = []
                for exon_id, exon_info in gene_hierarchy[gene_id]['mRNAs']\
                    [mRNA_id]['exons'].\
                    iteritems():

                    exon_rec = gene_hierarchy[gene_id]['mRNAs']\
                        [mRNA_id]['exons'][exon_id]['record']
                    strand = exon_rec.strand
                    exon_starts.append(exon_rec.start)
                    exon_ends.append(exon_rec.end)
                    mRNA.append(sorted([exon_rec.start, exon_rec.end]))

                mRNAs.append(mRNA)
            break

    mRNAs.sort(key=len)
    return tx_start, tx_end, exon_starts, exon_ends, gene_obj, \
           mRNAs, strand, chrom

