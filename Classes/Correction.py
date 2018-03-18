from Bio.Blast import NCBIXML
from copy import deepcopy

import logging

logger = logging.getLogger("root")

class Correction:
    def __init__(self, upstream_f, aligned_f, downstream_f):
        self.u = self.return_hsp(upstream_f)
        self.a = self.return_hsp(aligned_f)
        self.d = self.return_hsp(downstream_f)
        # This is the error value used to line up the upstream
        # and downstream sequences.
        # Initially setting to 7 (just over 2 AAs)
        self.error_distance = 7

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
        logger.debug("")
        if (self.u):
            logger.debug("ustream query range: %s:%s" % (self.u.query_start, self.u.query_end))
            logger.debug("ustream sbjct range: %s:%s" % (self.u.sbjct_start, self.u.sbjct_end))
        if (self.a):
            logger.debug("aligned query range: %s:%s" % (self.a.query_start, self.a.query_end))
            logger.debug("aligned sbjct range: %s:%s" % (self.a.sbjct_start, self.a.sbjct_end))
        if (self.d):
            logger.debug("dstream query range: %s:%s" % (self.d.query_start, self.d.query_end))
            logger.debug("dstream sbjct range: %s:%s" % (self.d.sbjct_start, self.d.sbjct_end))

    def process_logic_flag(self, flag, stats):
        logger.debug("Flag is: %s" % flag)
        if (flag == 1):
            stats.increment_upstream_aligned()
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
        self.debug_info()
        self.determine_relationships(stats)
