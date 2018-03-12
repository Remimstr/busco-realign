import os
from Bio import SeqIO
from Bio.Blast import NCBIXML

import logging

logger = logging.getLogger("root")

class Gene:
    def __init__(self, gene, path):
        base = os.path.basename(gene)
        self.name = os.path.splitext(base)[0]
        self.best_alignment = None
        self.record_files = self.split_records(gene, path)
        self.alignment_files = []

    def __str__(self):
        ret_str = "For gene %s: " % self.name
        if not self.best_alignment:
            ret_str += "There is no best record yet"
        else:
            ret_str += self.records[self.best].__str__()
        return ret_str

    def add_alignment_file(self, alignment_f):
        if (alignment_f not in self.alignment_files):
            self.alignment_files.append(alignment_f)

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
        return record_list

    # Chooses the best record that exists, remove it from the list of
    # record files and place it at "best"
    def choose_best_match(self):
        if not self.alignment_files:
            raise RuntimeError("Please make a list of records before choosing a best match")
        highest_b = 0
        best_alignment = None
        for alignment_f in self.alignment_files:
            with open(alignment_f, "rU") as result_h:
                try:
                    record = NCBIXML.read(result_h)
                    hsp = record.alignments[0].hsps[0]
                    b = hsp.bits
                    if (b > highest_b):
                        highest_b = b
                        best_alignment = alignment_f
                except:
                    logger.info("%s did not align, ignoring" % alignment_f)
        self.alignment_files.remove(best_alignment)
        self.best_alignment = best_alignment
        logger.info("%s has a new best alignment: %s" % (self.name, self.best_alignment))
