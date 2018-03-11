from Bio.Blast.Applications import NcbitblastnCommandline

import logging

logger = logging.getLogger("root")

def perform_alignment(subject, query, o_name, hsps=1):
    logger.info("Performing alignment on subject:\n\n%s\n\nand query: %s\n\nwith outfile: %s" % (subject, query, o_name))
    tblastn_cline = NcbitblastnCommandline(query=query, subject=subject, num_alignments=1, max_hsps=hsps, out=o_name, outfmt=5)
    tblastn_cline()
