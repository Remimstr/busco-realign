import os, sys

import logging

from alignment import *

logger = logging.getLogger("root")

# Manages the initial alignment where we align each gene with the query
def initial_alignment(container, directory):
    # Make a new directory for this process
    alignments_d = os.path.join(directory, "alignments")
    os.makedirs(alignments_d, exist_ok=True)

    for gene in container.genes:
        gene_alignments_d = os.path.join(alignments_d, gene.name)
        os.makedirs(gene_alignments_d, exist_ok=True)
        for key, value in gene.records.get_records().items():
            record_f = value.get_record_file_path()
            alignment_f = os.path.join(gene_alignments_d, key + ".xml")
            alignment_f = value.set_alignment_file_path(alignment_f)
            # And then process the new file
            if (not os.path.exists(alignment_f or container.overwrite)):
                perform_alignment(container.assembly.record, record_f, alignment_f)

        gene.set_best_record()
