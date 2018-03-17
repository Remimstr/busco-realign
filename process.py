import os, sys

import logging

from alignment import *

logger = logging.getLogger("root")

# Manages alignments and directory creation
def alignment(assembly, gene, stage_d, records):
    # Make a new directory for this process
    alignments_d = os.path.join(stage_d, "alignments")
    os.makedirs(alignments_d, exist_ok=True)

    gene_alignments_d = os.path.join(alignments_d, gene.name)
    os.makedirs(gene_alignments_d, exist_ok=True)
    for key, value in records.items():
        record_f = value.get_record_file_path()
        alignment_f = os.path.join(gene_alignments_d, key + ".xml")
        alignment_f = value.set_alignment_file_path(alignment_f)
        # And then process the new file
        if (not os.path.exists(alignment_f or container.overwrite)):
            perform_alignment(assembly, record_f, alignment_f)

# Split the best aligned records for further processing
def split_aligned_records(container, stage_two_d):
    # Make some new directories for this process
    alignments_d = os.path.join(stage_two_d, "alignments")
    os.makedirs(alignments_d, exist_ok=True)
    records_d = os.path.join(stage_two_d, "records")
    os.makedirs(records_d, exist_ok=True)

    for gene in container.genes:
        gene_records_d = os.path.join(records_d, gene.name)
        os.makedirs(gene_records_d, exist_ok=True)
        gene.best_record["data"].split_record(gene.fragments, gene_records_d)
