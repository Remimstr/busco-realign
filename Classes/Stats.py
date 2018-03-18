import logging

logger = logging.getLogger("root")

class Stats:
    def __init__(self):
        self.stats  = {
           "initial_alignments": {
               "both_aligned": 0,
               "upstream_only_aligned": 0,
               "downstream_only_aligned": 0,
               "neither_aligned": 0
           }
        }

    def increment_both_aligned(self):
        self.stats["initial_alignments"]["both_aligned"] += 1

    def increment_upstream_only_aligned(self):
        self.stats["initial_alignments"]["upstream_only_aligned"] += 1

    def increment_downstream_only_aligned(self):
        self.stats["initial_alignments"]["downstream_only_aligned"] += 1

    def increment_neither_aligned(self):
        self.stats["initial_alignments"]["neither_aligned"] += 1
