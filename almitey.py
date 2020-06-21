#! /usr/bin/env python3

import os
import sys
import argparse
import pickle

from kit.finder import SwitchpointFinder
from kit.utils import revcomp, remove_gaps_and_adjust_boundary
from kit.blast import get_blast_hits_with_orientation, quick_blastn
from kit.blast import set_blast_result_orientation, parse_blast_output

parser = argparse.ArgumentParser()

parser.add_argument(
    '-a', '--acefile', required=True
)

parser.add_argument(
    '-o', '--output-dir', required=False, default="./", type=str
)

parser.add_argument(
    '-w', '--window-size', required=False, default=7, type=int
)

parser.add_argument(
    '-d', '--min-depth', required=False, default=10, type=int
)

parser.add_argument(
    '-r', '--min-read-prop', required=False, default=0.01, type=float
)

args = parser.parse_args()

try:
    os.mkdir(args.output_dir)

except FileExistsError:
    pass

keyword_args = {
    'window_size': args.window_size,
    'min_depth': args.min_depth,
    'min_read_prop': args.min_read_prop
}

finder = SwitchpointFinder(args.acefile, args.output_dir, **keyword_args)
