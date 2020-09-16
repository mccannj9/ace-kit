#! /usr/bin/env python3

import argparse

from almitey import Almitey


parser = argparse.ArgumentParser()

parser.add_argument(
    '-i', '--input-dir', required=True
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

parser.add_argument(
    '-l', '--log-file', required=False, default=None
)

parser.add_argument(
    '-s', '--read-suffices', required=False, default="fr", type=str
)

args = parser.parse_args()

keyword_args = {
    'window_size': args.window_size,
    'min_depth': args.min_depth,
    'min_read_prop': args.min_read_prop
}

almitey_runner = Almitey(args.input_dir, args.output_dir, **keyword_args)

almitey_runner.run_on_all_clusters()
