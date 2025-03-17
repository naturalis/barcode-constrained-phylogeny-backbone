#!/usr/bin/env python
"""
Polytomy Resolution Tool - Main Script

A tool for resolving polytomies in phylogenetic trees using the Open Tree of Life (OpenToL) API.
This script serves as the command-line interface to the polytomy resolution pipeline.
"""

import sys
import argparse
import logging
import time
from polytomy.pipeline import PolytomyResolutionPipeline


# Set up logging
def setup_logging(log_level, log_file=None):
    """Configure logging system based on specified log level and optional log file."""
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level}")

    # Basic configuration for console logging
    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Add file handler if log_file is specified
    if log_file:
        # Create file handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(numeric_level)
        
        # Create formatter and add it to the handler
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        
        # Add file handler to the root logger
        logging.getLogger().addHandler(file_handler)
        
        logging.info(f"Logging to file: {log_file}")


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Resolve polytomies in phylogenetic trees using OpenToL API",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input tree file in Newick format"
    )

    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output resolved tree file"
    )

    parser.add_argument(
        "--alignment", "-a",
        help="Alignment file for branch length optimization"
    )

    parser.add_argument(
        "--sequences", "-s",
        help="FASTA file with additional sequences to place on the backbone"
    )

    parser.add_argument(
        "--optimization-tool",
        choices=["iqtree", "raxml-ng"],
        default="iqtree",
        help="Tool to use for branch length optimization"
    )

    parser.add_argument(
        "--max-memory",
        type=int,
        default=4000,
        help="Maximum memory usage in MB"
    )

    parser.add_argument(
        "--cache-dir",
        default=".opentol_cache",
        help="Directory to cache OpenToL API responses"
    )

    parser.add_argument(
        "--log-level",
        choices=["debug", "info", "warning", "error", "critical"],
        default="info",
        help="Set logging level"
    )

    parser.add_argument(
        "--skip-opentol",
        action="store_true",
        help="Skip OpenToL resolution and use only pruning-based resolution"
    )

    parser.add_argument(
        "--no-place-sequences",
        action="store_true",
        help="Skip placement of pruned and additional sequences"
    )

    parser.add_argument(
        "--version", "-v",
        action="version",
        version="%(prog)s 0.1.0"
    )

    parser.add_argument(
        "--match-tree-alignment",
        action="store_true",
        help="Match tree and alignment to ensure they have the same taxa"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=6,
        help="Number of threads for RAxML and IQTree"
    )

    parser.add_argument(
        "--model",
        default="GTRCAT",
        help="Sequence evolution model for RAxML"
    )

    parser.add_argument(
        "--compress-alignment",
        action="store_true",
        help="Compress alignment to most informative columns"
    )

    parser.add_argument(
        "--compress-columns",
        type=int,
        default=700,
        help="Number of columns to keep when compressing alignment"
    )

    parser.add_argument(
        "--filter-exemplars", 
        action="store_true",
        help="Filter tree to keep at most one exemplar per pair before optimization"
    )

    parser.add_argument(
        "--exemplar-pairs",
        help="Path to file with exemplar pair information (format: taxon<tab>exemplar1<tab>exemplar2)"
    )

    parser.add_argument(
        "--place-first-exemplars", 
        action="store_true",
        help="Place only first exemplar for each missing taxon"
    )

    parser.add_argument(
        "--graft-second-exemplars",
        action="store_true",
        help="Graft second exemplar of each pair next to the first one"
    )

    parser.add_argument(
        "--log-file",
        help="Path to output log file"
    )

    parser.add_argument(
        "--keep-files",
        action="store_true",
        help="Keep temporary files generated during sequence placement"
    )

    return parser.parse_args()


def main():
    """Main function."""
    # Parse command line arguments
    args = parse_args()

    # Set up logging
    setup_logging(args.log_level, args.log_file)
    logger = logging.getLogger(__name__)

    # Start time
    start_time = time.time()
    logger.info("Starting polytomy resolution")

    # Create configuration dict from arguments
    config = {
        'parser': {
            'max_memory': args.max_memory,
        },
        'opentol': {
            'cache_dir': args.cache_dir,
            'skip_opentol': args.skip_opentol,
        },
        'optimizer': {
            'tool': args.optimization_tool,
            'threads': args.threads,
            'model': args.model,
        },
        'placer': {
            'skip_placement': args.no_place_sequences,
            'threads': args.threads,
            'model': args.model,
            'compress': args.compress_alignment,
            'compress_columns': args.compress_columns,
            'keep_files': args.keep_files, 
        }
    }

    # Initialize pipeline
    pipeline = PolytomyResolutionPipeline(config=config)

    try:
        # Load input tree
        logger.info(f"Loading tree from {args.input}")
        pipeline.load_tree(args.input)

        # Match tree and alignment taxa if requested
        if args.match_tree_alignment and args.alignment:
            logger.info("Matching tree and alignment taxa")
            filtered_tree, filtered_alignment = pipeline.match_tree_and_alignment(
                args.input, 
                args.alignment
            )
    
            # Update paths if matching was successful
            if filtered_tree and filtered_alignment:
                args.input = filtered_tree
                args.alignment = filtered_alignment
                # Reload the tree with the filtered one
                pipeline.load_tree(filtered_tree)

        # Resolve polytomies (only if NOT using skip-opentol)
        if not args.skip_opentol:
            logger.info("Resolving polytomies")
            pipeline.resolve_polytomies()
        else:
            logger.info("Skipping polytomy resolution (--skip-opentol specified)")
            # Treat the input tree as already resolved
            pipeline.resolved_tree = pipeline.tree

            # If the tree is already optimized, set it directly
            if "optimized" in args.input:
                logger.info("Using pre-optimized input tree")
                pipeline.optimized_tree = pipeline.tree

        # Filter exemplar pairs if provided
        if args.exemplar_pairs and args.filter_exemplars:
            logger.info("Filtering to keep at most one exemplar per pair")
            pipeline.filter_exemplar_pairs(args.exemplar_pairs)

        # Optimize branch lengths if alignment provided and tree is not already optimized
        if args.alignment and not "optimized" in args.input:
            logger.info(f"Optimizing branch lengths using {args.optimization_tool}")
            pipeline.optimize_tree(args.alignment)
        else:
            logger.info("Skipping branch length optimization (tree already optimized)")

        # Place additional sequences if provided
        if args.sequences and not args.no_place_sequences:
            if args.exemplar_pairs and args.place_first_exemplars:
                logger.info(f"Placing first exemplars from {args.exemplar_pairs}")
                pipeline.place_additional_sequences(
                    args.sequences,
                    pairs_file=args.exemplar_pairs
                )
            else:
                logger.info(f"Placing all additional sequences from {args.sequences}")
                pipeline.place_additional_sequences(args.sequences)

        # Re-optimize branch lengths after placing first exemplars
        if args.exemplar_pairs and args.place_first_exemplars and args.alignment:
            logger.info("Re-optimizing branch lengths after placing first exemplars")
            pipeline.optimize_tree(args.alignment)

        # Graft second exemplars if provided
        if args.exemplar_pairs and args.graft_second_exemplars:
            logger.info(f"Grafting second exemplars from {args.exemplar_pairs}")
            pipeline.graft_second_exemplars(args.alignment, args.exemplar_pairs)

        # Re-optimize branch lengths after grafting second exemplars
        if args.exemplar_pairs and args.graft_second_exemplars and args.alignment:
            logger.info("Re-optimizing branch lengths after grafting")
            pipeline.reoptimize_tree_after_grafting(args.alignment)

        # Write output tree
        logger.info(f"Writing resolved tree to {args.output}")
        pipeline.write_tree(args.output)

        # Report statistics
        elapsed_time = time.time() - start_time
        logger.info(f"Polytomy resolution completed in {elapsed_time:.2f} seconds")

    except Exception as e:
        logger.error(f"Error during polytomy resolution: {str(e)}")
        logger.error("Exception details:", exc_info=True)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())