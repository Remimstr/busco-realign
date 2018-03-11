from Bio import SeqIO
import os

class Gene:
    def __init__(self, gene):
        self.name = os.path.basename(gene)
        self.records = gene
        self.best = None

    def __str__(self):
        ret_str = "For gene %s: " % self.name
        if (self.best == None):
            ret_str += "There is no best record yet"
        else:
            ret_str += self.records[self.best].__str__()
        return ret_str
