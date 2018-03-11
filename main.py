import argparse
import sys, os

import log

logger = log.setup_custom_logger("root")

from Classes import Container
from process import initial_alignment

description_message = """
A program for correcting read based on BUSCO alignments
"""

def process_args(args):
    gene_list = []
    if (args.buscoTSV):
        raise NotImplementedError("I haven't implemented this feature yet")
    elif (args.geneDirectory):
        genes = [os.path.join(args.geneDirectory, f) for f in os.listdir(args.geneDirectory)]
        genes = list(filter(lambda x: x.endswith(".fasta"), genes))
    logger.info("%s genes will be processed" % len(genes))

    # Prepares a new directory for intermediate files
    directory = "tmp"
    os.makedirs(directory, exist_ok=True)

    container = Container(args.assembly, genes, os.path.join(os.getcwd(), directory), args.forceOverwrite)

    initial_alignment(container, directory)

def main():
    parser = argparse.ArgumentParser(description=description_message)
    parser.add_argument("-a", "--assembly", action="store", required=True, type=str)
    mxg = parser.add_mutually_exclusive_group(required=True)
    mxg.add_argument("-d", "--geneDirectory", action="store", type=str)
    mxg.add_argument("-t", "--buscoTSV", action="store", type=str)
    parser.add_argument("-o", "--forceOverwrite", action="store", type=bool, default=False)
    args = parser.parse_args()
    args.assembly = "/mnt/extra_storage/jshlorv/minION_data/psy_508_seqRun/assemblies/miniasm_test_assemblies/508_v2/qu75/508_v2_pc_mLmH_GtQu75_tr50.1000_1000_02_50.contigs.ONT.fa"
    args.geneDirectory = "../508_v2_pc_mLmH_GtQu75_tr50.1000_1000_02_50.contigs.ONT_gam_1e-03/Fragmented"
    process_args(args)

if __name__ == "__main__":
    main()
