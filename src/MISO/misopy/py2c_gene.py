import misopy
import pysplicing

def py2c_gene(py_gene):
    """
    Convert a Python Gene object to a C gene object for use
    with C MISO.
    """
    # Description of exon lens
    CMISO_exon_lens = tuple([(part.start, part.end) \
                             for part in py_gene.parts])

    # Description of isoforms of gene
    isoforms_desc = []
    for isoform in py_gene.isoforms:
        curr_iso_desc = tuple([py_gene.parts.index(iso_part) \
                               for iso_part in isoform.parts])
        isoforms_desc.append(curr_iso_desc)

    CMISO_isoforms_desc = tuple(isoforms_desc)
    c_gene = pysplicing.createGene(CMISO_exon_lens,
                                   CMISO_isoforms_desc)
    return c_gene
