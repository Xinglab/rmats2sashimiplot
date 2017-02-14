##
## Test Pysplicing (fastmiso)
##

import misopy
import misopy.gff_utils
import pysplicing

gene = pysplicing.createGene( ((1,100), (201,300), (401,500)),
                              ((0,1), (0,2), (0,1,2)) )

pysplicing.noIso(gene)
pysplicing.isoLength(gene)

reads = pysplicing.simulateReads(gene, 0L, (0.2,0.3,0.5), 2000L, 33L)

# Load BAM file


print reads[1], type(reads[1])
results = pysplicing.MISO(gene, 0L, reads[1], reads[2], 33L, 5000L, 500L, 10L,
                          (1.0, 1.0, 1.0))


# est2 = pysplicing.solveIsoGene(gene, 0L, 33L, reads[1], reads[2])
########################
