from Bio.Blast import NCBIXML
from copy import deepcopy

import logging

logger = logging.getLogger("root")

class Correction:
    def __init__(self, upstream_f, aligned_f, downstream_f):
        self.upstream_o = self.return_hsp(upstream_f)
        self.aligned_o = self.return_hsp(aligned_f)
        self.downstream_o = self.return_hsp(downstream_f)

    # This function opens an alignment and returns the hsp object
    def return_hsp(self, alignment_f):
        with open(alignment_f, "rU") as alignment_h:
            try:
                record = deepcopy(NCBIXML.read(alignment_h))
                hsp = record.alignments[0].hsps[0]
                return hsp
            except:
                return None

    def correct(self, container):
        if not self.aligned_o:
            logger.error("An re-alignment of an aligned region was not found, this is bad")
        elif (self.upstream_o and self.downstream_o):
            container.stats.increment_both_aligned()
        elif (self.upstream_o):
            container.stats.increment_upstream_only_aligned()
        elif (self.downstream_o):
            container.stats.increment_downstream_only_aligned()
        else:
            container.stats.increment_neither_aligned()
