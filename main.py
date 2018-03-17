import argparse
import sys, os
import pickle

import log

logger = log.setup_custom_logger("root")
pickle_file = "container_pickle.pkl"

from Classes import Container
from process import initial_alignment

description_message = """
A program for correcting read based on BUSCO alignments
"""

# Use pickle to persist classes across stages (for easier development)
def load_pickle(pickle_path):
    with open(pickle_path, "rb") as pp:
        return pickle.load(pp)

def save_pickle(container, pickle_path):
    with open(pickle_path, "wb") as pp:
        pickle.dump(container, pp, pickle.HIGHEST_PROTOCOL)

def stage_one(container, directory):
    stage_one_d = os.path.join(directory, "stage_one")
    os.makedirs(stage_one_d, exist_ok=True)
    initial_alignment(container, stage_one_d)

def process_args(args):
    gene_list = []
    if (args.buscoTSV):
        raise NotImplementedError("I haven't implemented this feature yet")
    elif (args.geneDirectory):
        genes = [os.path.join(args.geneDirectory, f) for f in os.listdir(args.geneDirectory)]
        genes = list(filter(lambda x: x.endswith(".fasta"), genes))
    logger.info("%s genes will be processed" % len(genes))
    return genes

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
    genes = process_args(args)

    # Prepares a new directory for intermediate files
    tmp = "tmp"
    os.makedirs(tmp, exist_ok=True)

    # Make a new pickle path for class storage
    pickle_path = os.path.join(tmp, pickle_file)

    if not os.path.exists(pickle_path):
        container = Container(args.assembly, genes, os.path.join(os.getcwd(), tmp), args.forceOverwrite)
        stage_one(container, tmp)
        save_pickle(container, pickle_path)
    else:
        container = load_pickle(pickle_path)


if __name__ == "__main__":
    main()
