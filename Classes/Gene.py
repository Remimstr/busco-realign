from Bio import SeqIO

class Gene:
    def __init__(self, gene):
        gene_obj = list(SeqIO.parse(gene, "fasta"))

    def __str__(self):
        print "A gene will be here"
