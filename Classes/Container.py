import os
from copy import deepcopy

from .Assembly import Assembly
from .Gene import Gene

class Container:
    def __init__(self, assembly, out_assembly, genes, path, overwrite):
        #print(Gene(genes[0]))
        self.assembly = Assembly(assembly)
        self.out_assembly = out_assembly
        self.genes = [Gene(gene, path) for gene in genes]
        self.overwrite = overwrite
        self.savefile = os.path.join(path, "dump.txt")
        self.corrections = []

    def __str__(self):
        return """
        Assembly: %s\n
        There are %s genes
        """ % (self.assembly, len(self.genes))

    def add_correction(self, correction):
        self.corrections.append(correction)

    def copy_assembly(self):
        return deepcopy(self.assembly.record)
