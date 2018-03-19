from Bio.Blast import NCBIXML
from Bio import SeqIO
from copy import deepcopy

import logging

logger = logging.getLogger("root")

class Correction:
    def __init__(self, assembly_f, upstream_f, aligned_f, downstream_f, full_gene_f):
        self.assembly = assembly_f
        self.u = self.return_hsp(upstream_f)
        self.a = self.return_hsp(aligned_f)
        self.d = self.return_hsp(downstream_f)
        self.f = self.return_hsp(full_gene_f)
        # This is the error value used to line up the upstream
        # and downstream sequences.
        # Initially setting to 7 (just over 2 AAs)
        self.error_distance = 7

    # Returns the sequence of assembly between start and end
    def index_assembly(self, assembly, start, end):
        with open(assembly, "rU") as assembly_h:
            assembly_o = SeqIO.read(assembly_h, "fasta")
            return assembly_o.seq[start:end].translate()

    def return_assembly_copy(self):
        return deepcopy(self.assembly)

    # This function opens an alignment and returns the hsp object
    def return_hsp(self, alignment_f):
        with open(alignment_f, "rU") as alignment_h:
            try:
                record = deepcopy(NCBIXML.read(alignment_h))
                hsp = record.alignments[0].hsps[0]
                return hsp
            except:
                return None

    def collect_initial_stats(self, stats):
        if not self.a:
            logger.error("An re-alignment of an aligned region was not found, this is bad")
            stats.increment_errors()
        elif (self.u and self.d):
            stats.increment_both_aligned()
        elif (self.u):
            stats.increment_upstream_only_aligned()
        elif (self.d):
            stats.increment_downstream_only_aligned()
        else:
            stats.increment_neither_aligned()

    def debug_info(self):
        if (self.u):
            logger.debug("ustream query range: %s:%s" % (self.u.query_start, self.u.query_end))
            logger.debug("ustream sbjct range: %s:%s" % (self.u.sbjct_start, self.u.sbjct_end))
        if (self.a):
            logger.debug("aligned query range: %s:%s" % (self.a.query_start, self.a.query_end))
            logger.debug("aligned sbjct range: %s:%s" % (self.a.sbjct_start, self.a.sbjct_end))
        if (self.d):
            logger.debug("dstream query range: %s:%s" % (self.d.query_start, self.d.query_end))
            logger.debug("dstream sbjct range: %s:%s" % (self.d.sbjct_start, self.d.sbjct_end))


    ### REPAIR SECTION

    # This is where the magic happens. Take two positions and try and resolve
    # the base pairs on the assembly between them.

    def determine_num_mutations(self, pos1, pos2):
        # This is some weird logic to circumvent return values of 0 or 3 (would do nothing)
        offset = ((abs(pos1 - pos2) + 1) % 3)
        if (offset != 0) and (offset != 3):
            return offset
        else:
            return 2

    # Provides useful logging messages for determining if the correction is working
    def exploratory(self, first_a, second_a, frame_offset):
        first_query_seq = first_a.query[:first_a.query_end]
        second_query_seq = second_a.query[second_a.query_start:]
        first_sbjct_start = first_a.sbjct_start
        first_sbjct_end = first_a.sbjct_end
        second_sbjct_start = second_a.sbjct_start
        second_sbjct_end = second_a.sbjct_end

        logger.debug("Distance          : %s" % abs(first_sbjct_end - second_sbjct_start))
        logger.debug("Frame Offset      : %s" % frame_offset)
        logger.debug("First Query Seq   : %s" % first_query_seq)
        logger.debug("Second Query Seq  : %s" % second_query_seq)
        logger.debug("First Query Frame : %s" % self.index_assembly(self.assembly, first_sbjct_start-1, first_sbjct_end))
        logger.debug("Second Query Frame: %s" % self.index_assembly(self.assembly, second_sbjct_start-1, second_sbjct_end))
        logger.debug("Frame 1           : %s" % self.index_assembly(self.assembly, first_sbjct_start-1, second_sbjct_end))
        logger.debug("Frame 2           : %s" % self.index_assembly(self.assembly, first_sbjct_start-2, second_sbjct_end))
        logger.debug("Frame 3           : %s" % self.index_assembly(self.assembly, first_sbjct_start-3, second_sbjct_end))

    def correct_frame(self, first_a, second_a):
        # It seems counter-intuitive but we want to fix the end of the first alignment and the start of the second
        first_sbjct_start = first_a.sbjct_start
        first_sbjct_end = first_a.sbjct_end
        second_sbjct_start = second_a.sbjct_start
        second_sbjct_end = second_a.sbjct_end
        frame_offset = self.determine_num_mutations(first_sbjct_end, second_sbjct_start)
        self.exploratory(first_a, second_a, frame_offset)

        logger.debug("*** BEGIN CORRECTION OF FRAME %s ***" % frame_offset)
        # Get a new copy of the assembly
        with open (self.return_assembly_copy(), "rU") as assembly_h:
            assembly = SeqIO.read(assembly_h, "fasta")
            # Add a new base pair
            assembly.seq = assembly.seq[:first_sbjct_end-1] + (frame_offset * "a") + assembly.seq[second_sbjct_start-1:]
            logger.debug("Added %s new base pairs" % frame_offset)
            logger.debug("New Frame         : %s" % (assembly.seq[first_sbjct_start-1:second_sbjct_end]).translate())
            logger.debug("Old Alignment     : %s" % (self.f.query))
        logger.debug("*** END CORRECTION OF FRAME %s ***" % frame_offset)

    ### END OF REPAIR SECTION

    def process_logic_flag(self, flag, stats):
        # logger.debug("Flag is: %s" % flag)
        if (flag == 1):
            stats.increment_upstream_aligned()
            self.correct_frame(self.u, self.a)

        elif (flag == 3):
            stats.increment_aligned_downstream()
        elif (flag == 4):
            stats.increment_upstream_aligned_downstream()
        elif (flag == 5):
            stats.increment_aligned_upstream()
        elif (flag == 7):
            stats.increment_downstream_aligned()
        elif (flag == 12):
            stats.increment_downstream_aligned_upstream()
        elif (flag == 0):
            stats.increment_no_relationship()
        else:
            logger.debug("Undetermined flag: %s" % flag)
            stats.increment_undetermined_relationship()

    def determine_relationships(self, stats):
        # Use some flag logic to resolve all states
        flag = 0
        if (self.u and self.a):
            # Case where upstream end is upstream of the alignment start
            if (abs(self.a.sbjct_start - self.u.sbjct_end) < self.error_distance):
                flag += 1
            # Case where the upstream start is dowmstream of the alignment end
            if (abs(self.a.sbjct_end - self.u.sbjct_start) < self.error_distance):
                flag += 5
        if (self.a and self.d):
            # Case where downstream start is downstream of the alignment end
            if (abs(self.d.sbjct_start - self.a.sbjct_end) < self.error_distance):
                flag += 3
            # Case where the downstream end is upstream of the alignment start
            if (abs(self.d.sbjct_end - self.a.sbjct_start) < self.error_distance):
                flag += 7
        self.process_logic_flag(flag, stats)

    def correct(self, stats):
        self.collect_initial_stats(stats)
        # self.debug_info()
        self.determine_relationships(stats)
