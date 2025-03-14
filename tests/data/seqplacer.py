from polytomy.sequence_placer import SequencePlacer
from polytomy.tree_parser import TreeParser
import argparse

# parse command line arguments: --tree, --alignment, --output, --threads, --model, --prefix, --keep-files, --output-dir
parser = argparse.ArgumentParser()
parser.add_argument("--intree", "-i", required=True)
parser.add_argument("--alignment", "-a", required=True)
parser.add_argument("--output", "-o", required=True)
parser.add_argument("--threads", "-t", type=int, default=6)
parser.add_argument("--model", "-m", default="GTRCAT")
parser.add_argument("--prefix", "-p", default="seq_placement")
parser.add_argument("--keep-files", "-k", action="store_true")
parser.add_argument("--dir", "-d", default=".")
args = parser.parse_args()

parser = TreeParser()
backbone_tree = parser.parse_from_file(args.intree)
placer = SequencePlacer(backbone_tree, config={
    'threads': args.threads,
    'model': args.model,
    'prefix': args.prefix,
    'keep_files': args.keep_files,
    'output_dir': args.dir
})
result_tree = placer.place_sequences(args.alignment, prefilter=True, compress=True)