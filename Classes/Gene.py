import os
from Bio import SeqIO

import logging
from .Records import Records

logger = logging.getLogger("root")

class Gene:
    def __init__(self, gene, path):
        base = os.path.basename(gene)
        self.name = os.path.splitext(base)[0]
        self.best_record = {"name": None, "data": None}
        self.fragments = {"upstream": None, "aligned": None, "downstream": None}
        self.records = Records()
        self.split_records(gene, path)

    def __str__(self):
        ret_str = "For gene %s: " % self.name
        if not self.best_alignment:
            ret_str += "There is no best record yet"
        else:
            ret_str += self.records[self.best].__str__()
        return ret_str

    # Takes a path, and a gene object and makes a new directory
    # for the gene at path, then splits the gene's records and places
    # them inside that directory
    def split_records(self, gene, path):
        directory = os.path.join(path, "records")
        os.makedirs(directory, exist_ok=True)
        record_list = []
        gene_d = os.path.join(directory, self.name)
        os.makedirs(gene_d, exist_ok=True)
        for record in list(SeqIO.parse(gene, "fasta")):
            record_fasta_f = os.path.join(gene_d, record.name + ".fasta").replace(":", "-")
            # We need to create a new file for each record we just found
            if (not os.path.exists(record_fasta_f)):
                with open(record_fasta_f, "w") as output_h:
                    SeqIO.write(record, output_h, "fasta")
            self.records.create_record(record.name, record_fasta_f)
        logger.info("Split gene %s into %s records" % (self.name, len(self.records.get_records())))

    def set_best_record(self):
        self.best_record = self.records.get_best_record(self.best_record)
        logger.info("Chose a new best record %s for gene %s" % (self.best_record["name"], self.name))
