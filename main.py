import argparse
import logging, sys

import os

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

description_message = """
A program for correcting read based on BUSCO alignments
"""

def process_args(args):
    gene_list = []
    if (args.buscoTSV):
        raise NotImplementedError("I haven't implemented this feature yet")
    elif (args.geneDirectory):
        gene_list = filter(lambda x: x.endswith(".fasta"), os.listdir(args.geneDirectory))
    logging.info("%s genes will be processed" % len(gene_list))

def main():
    parser = argparse.ArgumentParser(description=description_message)
    parser.add_argument("-a", "--assembly", action="store", required=True, type=str)
    mxg = parser.add_mutually_exclusive_group(required=True)
    mxg.add_argument("-d", "--geneDirectory", action="store", type=str)
    mxg.add_argument("-t", "--buscoTSV", action="store", type=str)
    args = parser.parse_args()
    args.assembly = "/mnt/extra_storage/jshlorv/minION_data/psy_508_seqRun/assemblies/miniasm_test_assemblies/508_v2/qu75/508_v2_pc_mLmH_GtQu75_tr50.1000_1000_02_50.contigs.ONT.fa"
    args.geneDirectory = "../508_v2_pc_mLmH_GtQu75_tr50.1000_1000_02_50.contigs.ONT_gam_1e-03/Fragmented"
    process_args(args)

if __name__ == "__main__":
    main()
