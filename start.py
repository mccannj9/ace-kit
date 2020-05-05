#! /usr/bin/env python3

import os
import argparse

from matplotlib import pyplot

from kit.finder import SwitchpointFinder

parser = argparse.ArgumentParser()

parser.add_argument('-a', '--acefile', required=True)
parser.add_argument('-o', '--output-dir', required=False, default="./", type=str)
parser.add_argument('-w', '--window-size', required=False, default=7, type=int)
parser.add_argument('-d', '--min-depth', required=False, default=10, type=int)
parser.add_argument('-r', '--min-read-prop', required=False, default=0.1, type=float)
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
acefile = finder.acefile

results_dict = finder.fit()
