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
def setup_logging(log_level):
    """Configure logging system based on specified log level."""
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level}")

    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


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

    return parser.parse_args()


def main():
    """Main function."""
    # Parse command line arguments
    args = parse_args()

    # Set up logging
    setup_logging(args.log_level)
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
        },
        'placer': {
            'skip_placement': args.no_place_sequences,
        }
    }

    # Initialize pipeline
    pipeline = PolytomyResolutionPipeline(config=config)

    try:
        # Load input tree
        logger.info(f"Loading tree from {args.input}")
        pipeline.load_tree(args.input)

        # Resolve polytomies
        logger.info("Resolving polytomies")
        pipeline.resolve_polytomies()

        # Optimize branch lengths if alignment provided
        if args.alignment:
            logger.info(f"Optimizing branch lengths using {args.optimization_tool}")
            pipeline.optimize_tree(args.alignment)

        # Place additional sequences if provided
        if args.sequences and not args.no_place_sequences:
            logger.info(f"Placing additional sequences from {args.sequences}")
            pipeline.place_additional_sequences(args.sequences)

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