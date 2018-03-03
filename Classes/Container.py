from Assembly import Assembly
from Gene import Gene

class Container:
    def __init__(self, assembly, genes):
        #print(Gene(genes[0]))
        self.assembly = Assembly(assembly)
        self.genes = [Gene(gene) for gene in genes]

    def __str__(self):
        return """
        Assembly: %s\n
        There are %s genes
        """ % (self.assembly, len(self.genes))
