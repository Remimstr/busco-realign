from Bio import SeqIO

class Assembly:
    def __init__ (self, assembly):
        self.record = assembly

    def __str__(self):
        return self.record
