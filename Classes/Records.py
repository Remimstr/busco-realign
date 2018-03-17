import os
from copy import deepcopy
from Bio import SeqIO
from Bio.Blast import NCBIXML
import logging

logger = logging.getLogger("root")

class Record:
    def __init__(self, record_file_path="", alignment_file_path=""):
        self.data = {"record": record_file_path, "alignment": alignment_file_path}

    def __str__(self):
        return self.data

    def get_record_file_path(self):
        return self.data["record"]

    def get_alignment_file_path(self):
        return self.data["alignment"]

    def set_record_file_path(self, record_file_path):
        self.data["record"] = record_file_path.replace(":", "-")
        return self.data["record"]

    def set_alignment_file_path(self, alignment_file_path):
        self.data["alignment"] = alignment_file_path.replace(":", "-")
        return self.data["alignment"]

    def write_to_file(self, record_o, record_f):
        with open(record_f, "w") as record_h:
           SeqIO.write(record_o, record_h, "fasta")
        return record_f

    # Returns a subset of a record based on positional information
    def subset_record(self, record_f, start, end=None):
        record = None
        with open(record_f, "rU") as record_h:
            record = deepcopy(SeqIO.read(record_h, "fasta"))
            if not end:
                record.seq = record.seq[start:]
            else:
                record.seq = record.seq[start:end]
        return record

    # Splits a record based on alignment
    def split_record(self, fragments, path):
        with open(self.get_alignment_file_path(), "rU") as alignment_h:
            record = NCBIXML.read(alignment_h).alignments[0].hsps[0]
            record_f = self.get_record_file_path()
            # Process upstream
            upstream_o = self.subset_record(record_f, 0, record.query_start)
            upstream_f = os.path.join(path, "upstream.fasta")
            self.write_to_file(upstream_o, upstream_f)
            fragments["upstream"] = Record(upstream_f)
            # Process aligned
            aligned_o = self.subset_record(record_f, record.query_start, record.query_end)
            aligned_f = os.path.join(path, "aligned.fasta")
            self.write_to_file(aligned_o, aligned_f)
            fragments["aligned"] = Record(aligned_f)
            # Process downstream
            downstream_o = self.subset_record(record_f, record.query_end)
            downstream_f = os.path.join(path, "downstream.fasta")
            self.write_to_file(downstream_o, downstream_f)
            fragments["downstream"] = Record(downstream_f)

class Records:
    def __init__(self):
        self.records = {}

    def check_records(self, op, value = None):
        if (not value and not self.records):
            raise RuntimeError("Please make a list of records before %s" % op)
        elif value:
            try:
                self.records[value]
            except:
                raise RuntimeError("Object %s is not found in record %s" % value)

    def get_records(self):
        return self.records

    def get_record(self, index):
        return self.records[index]

    def get_record_paths(self):
        self.check_records("getting record paths")
        records = []
        for (key, value) in self.get_records():
            records.append(value["record"])
        return records

    def get_alignment_paths(self):
        self.check_records("getting alignment paths")
        alignments = []
        for (key, value) in self.get_records():
            alignments.append(value["alignment"])
        return alignments

    def create_record(self, name, record_file_path, alignment_file_path):
        new_record = Record(record_file_path, alignment_file_path)
        self.records[name] = new_record

    # Chooses the best record that exists, remove it from the list of
    # record files and place it at "best"
    def get_best_record(self, best_record):
        self.check_records("choosing a best match")
        highest_b = 0
        for key, value in self.get_records().items():
            alignment_f = value.get_alignment_file_path()
            with open(alignment_f, "rU") as result_h:
                try:
                    record = NCBIXML.read(result_h)
                    hsp = record.alignments[0].hsps[0]
                    b = hsp.bits
                    if (b > highest_b):
                        highest_b = b
                        best_record["name"] = key
                        best_record["data"] = value
                except:
                    logger.info("%s did not align, ignoring" % key)
        del self.records[best_record["name"]]
        return best_record
