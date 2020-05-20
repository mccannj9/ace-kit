#! /usr/bin/env python3

import os
import argparse

from kit.finder import SwitchpointFinder
from kit.utils import revcomp

parser = argparse.ArgumentParser()

parser.add_argument('-a', '--acefile', required=True)

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
acefile = finder.acefile

results_dict, all_reads = finder.fit()

results = finder.results

# boundary_reads = [results[k].reads for k in results]
# boundary_reads = [item for sublist in boundary_reads for item in sublist]

boundary_read_prefixes = set([y.name[:-1] for x,y in all_reads.items() if y.boundary])

mate_pairs = []
for prefix in boundary_read_prefixes:
    f = f"{prefix}f"
    r = f"{prefix}r"
    if (f in all_reads) and (r in all_reads):
        rf = all_reads[f]
        rr = all_reads[r]
        mate_pairs.append(
            (rf.boundary * rr.boundary, rf, rr)
        )

boundary_mate_pairs = [(x,y) for w,x,y in mate_pairs if w]

for x, y in boundary_mate_pairs:
    print(f">{x.name}\n{x.seq.replace('*', '-')}")
    print(f">{y.name}\n{revcomp(y.seq)}")
