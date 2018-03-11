import os, sys
from Bio.Blast import NCBIXML
from Bio import SeqIO

import logging

from alignment import *

logger = logging.getLogger("root")

# Chooses the best match from a series of alignments and records it
# in the gene object
def choose_best_match(gene_obj, xml):
    highest_b = 0
    best_record = None
    with open(xml, "rU") as result_h:
        for record in NCBIXML.parse(result_h):
            # Let's assume only hsp per iteration
            hsp = record.alignments[0].hsps[0]
            print(record.alignments[0].hit_id)
            b = hsp.bits
            if (b > highest_b):
                highest_b = b
    sys.exit(0)


# Manages the initial alignment where we align each gene with the query
def initial_alignment(container, directory):
    # Make a new directory for this process
    """
    records_d = os.path.join(directory, "records")
    alignments_d = os.path.join(directory, "alignments")
    os.makedirs(records_d, exist_ok=True)
    os.makedirs(alignments_d, exist_ok=True)

    for gene in container.genes:
        gene_records_d = os.path.join(records_d, gene.name)
        gene_alignments_d = os.path.join(alignments_d, gene.name)
        os.makedirs(gene_records_d, exist_ok=True)
        os.makedirs(gene_alignments_d, exist_ok=True)
        for record in list(SeqIO.parse(gene.records, "fasta")):
            alignment_fasta_f = os.path.join(gene_alignments_d, record.name + ".xml")
            # And then process the new file
            if (not os.path.exists(alignment_fasta_f) or container.overwrite):
                perform_alignment(container.assembly.record, record_fasta_f, alignment_fasta_f)
        #choose_best_match(gene, o_file)
    """
