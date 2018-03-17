from Bio.Blast import NCBIXML
import logging

logger = logging.getLogger("root")

class Record:
    def __init__(self, record_file_path, alignment_file_path):
        self.data = {"record": record_file_path, "alignment": alignment_file_path}

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

    def create_record(self, name, record_file_path = "", alignment_file_path = ""):
        new_record = Record(record_file_path, alignment_file_path)
        self.records[name] = new_record

    # Chooses the best record that exists, remove it from the list of
    # record files and place it at "best"
    def get_best_record(self):
        self.check_records("choosing a best match")
        highest_b = 0
        best_record = None
        for key, value in self.get_records().items():
            alignment_f = value.get_alignment_file_path()
            with open(alignment_f, "rU") as result_h:
                try:
                    record = NCBIXML.read(result_h)
                    hsp = record.alignments[0].hsps[0]
                    b = hsp.bits
                    if (b > highest_b):
                        highest_b = b
                        best_record = key
                except:
                    logger.info("%s did not align, ignoring" % key)
        del self.records[best_record]
        return best_record
