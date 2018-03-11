import os
from Bio import SeqIO

import logging

logger = logging.getLogger("root")

class Gene:
    def __init__(self, gene, path):
        base = os.path.basename(gene)
        self.name = os.path.splitext(base)[0]
        self.best = None
        self.record_files = self.split_records(gene, path)

    def __str__(self):
        ret_str = "For gene %s: " % self.name
        if (self.best == None):
            ret_str += "There is no best record yet"
        else:
            ret_str += self.records[self.best].__str__()
        return ret_str

    # Takes a path, and a gene object and makes a new directory
    # for the gene at path, then splits the gene's records and places
    # them inside that directory
    def split_records(self, gene, path):
        record_list = []
        gene_d = os.path.join(path, self.name)
        os.makedirs(gene_d, exist_ok=True)
        for record in list(SeqIO.parse(gene, "fasta")):
            record_fasta_f = os.path.join(gene_d, record.name + ".fasta").replace(":", "-")
            # We need to create a new file for each record we just found
            if (not os.path.exists(record_fasta_f)):
                with open(record_fasta_f, "w") as output_h:
                    SeqIO.write(record, output_h, "fasta")
            record_list.append(record_fasta_f)
        logger.info("Split gene %s into %s records" % (self.name, len(record_list)))
