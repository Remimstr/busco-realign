import json

import logging

logger = logging.getLogger("root")

class Stats:
    def __init__(self):
        self.window_size = 5
        self.kmer_size = 4
        # Dictionary of kmer frequencies
        self.kmers = {}
        self.stats  = {
           "initial_alignments": {
               "both_aligned": 0,
               "upstream_only_aligned": 0,
               "downstream_only_aligned": 0,
               "neither_aligned": 0
           },
            "relationships": {
                "none": 0,
                "undetermined": 0,
                "forward": {
                    "upstream-aligned-downstream": 0,
                    "upstream-aligned": 0,
                    "aligned-downstream": 0
                },
                "reverse": {
                    "downstream-aligned-upstream": 0,
                    "downstream-aligned": 0,
                    "aligned-upstream": 0
                },
            },
            "errors": 0
        }

    def __str__(self):
        ret = "\nSummary Stats: \n"
        ret += json.dumps(self.stats, sort_keys=True, indent=4)
        return ret

    def increment_errors(self):
        self.stats["errors"] += 1

    def increment_both_aligned(self):
        self.stats["initial_alignments"]["both_aligned"] += 1

    def increment_upstream_only_aligned(self):
        self.stats["initial_alignments"]["upstream_only_aligned"] += 1

    def increment_downstream_only_aligned(self):
        self.stats["initial_alignments"]["downstream_only_aligned"] += 1

    def increment_neither_aligned(self):
        self.stats["initial_alignments"]["neither_aligned"] += 1

    def increment_undetermined_relationship(self):
        self.stats["relationships"]["undetermined"] += 1

    def increment_no_relationship(self):
        self.stats["relationships"]["none"] += 1

    def increment_upstream_aligned_downstream(self):
        self.stats["relationships"]["forward"]["upstream-aligned-downstream"] += 1

    def increment_upstream_aligned(self):
        self.stats["relationships"]["forward"]["upstream-aligned"] += 1

    def increment_aligned_downstream(self):
        self.stats["relationships"]["forward"]["aligned-downstream"] += 1

    def increment_downstream_aligned_upstream(self):
        self.stats["relationships"]["reverse"]["downstream-aligned-upstream"] += 1

    def increment_downstream_aligned(self):
        self.stats["relationships"]["reverse"]["downstream-aligned"] += 1

    def increment_aligned_upstream(self):
        self.stats["relationships"]["reverse"]["aligned-upstream"] += 1

    ### KMER STATS SECTION
    def add_kmer(self, kmer):
        if kmer in self.kmers:
            self.kmers[kmer] += 1
        else:
            self.kmers[kmer] = 1

    def dump_kmer_stats(self, kmer_stats_f):
        with open(kmer_stats_f, "w") as outfile:
            json.dump(self.kmers, outfile)
        logger.info("Dumped kmer stats to file %s" % kmer_stats_f)
