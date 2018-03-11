import os

from .Assembly import Assembly
from .Gene import Gene

class Container:
    def __init__(self, assembly, genes, path, overwrite):
        #print(Gene(genes[0]))
        self.assembly = Assembly(assembly)
        directory = os.path.join(path, "records")
        os.makedirs(directory, exist_ok=True)
        self.genes = [Gene(gene, directory) for gene in genes]
        self.overwrite = overwrite

    def __str__(self):
        return """
        Assembly: %s\n
        There are %s genes
        """ % (self.assembly, len(self.genes))
