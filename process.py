import os, sys
from Bio import SeqIO

import logging

from alignment import perform_alignment

from Classes import Correction
from Classes import Stats

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

# Manages the correction stage
def correction(container, stage_three_d):
    stats = Stats()
    assembly_f = container.copy_assembly()
    with open(container.copy_assembly(), "r") as assembly_f:
        assembly_o = SeqIO.read(assembly_f, "fasta")
        for gene in container.genes:
            upstream = gene.fragments["upstream"].get_alignment_file_path()
            aligned = gene.fragments["aligned"].get_alignment_file_path()
            downstream = gene.fragments["downstream"].get_alignment_file_path()
            correction = Correction(container.assembly.record, upstream, aligned, downstream, gene.best_record["data"].get_alignment_file_path())
            correction.correct(assembly_o, stats)
            container.add_correction(correction)
        stats.dump_kmer_stats("kmer_stats.json")
        print(stats)
    with open(container.out_assembly, "w") as assembly_h:
        SeqIO.write(assembly_o, assembly_h, "fasta")
