from Bio import SeqIO

class Assembly:
    def __init__ (self, assembly):
        self.assembly = SeqIO.read(assembly, "fasta")

    def __str__(self):
        return "An assembly will be here"
