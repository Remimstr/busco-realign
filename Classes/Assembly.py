from Bio import SeqIO

class Assembly:
    def __init__ (self, assembly):
        self.record = SeqIO.read(assembly, "fasta")

    def __str__(self):
        return self.record.__str__()
