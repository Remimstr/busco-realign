import os

from .Assembly import Assembly
from .Gene import Gene
from .Stats import Stats

class Container:
    def __init__(self, assembly, genes, path, overwrite):
        #print(Gene(genes[0]))
        self.assembly = Assembly(assembly)
        self.genes = [Gene(gene, path) for gene in genes]
        self.overwrite = overwrite
        self.savefile = os.path.join(path, "dump.txt")
        self.corrections = []
        self.stats = Stats()

    def __str__(self):
        return """
        Assembly: %s\n
        There are %s genes
        """ % (self.assembly, len(self.genes))

    def add_correction(self, correction):
        self.corrections.append(correction)
