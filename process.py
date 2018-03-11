import os

import logging

from alignment import *

logger = logging.getLogger("root")

def initial_alignment(container, directory):
    # Make a new directory for this process
    subdir = "stage_one"
    subdir = os.path.join(directory, subdir)
    os.makedirs(subdir, exist_ok=True)

    for gene in container.genes:
        o_file = os.path.join(subdir, gene.name)
        perform_alignment(container.assembly.record, gene.records, o_file)
