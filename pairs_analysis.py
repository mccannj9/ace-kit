#! /usr/bin/env python3

import subprocess
import argparse
import sys
import pickle

from difflib import get_close_matches

from kit.contig import Read, Pair


parser = argparse.ArgumentParser()

parser.add_argument('-p', '--pickle', required=True, type=str)

parser.add_argument(
    '-o', '--output-dir', required=False, default="./", type=str
)

parser.add_argument(
    '-k', '--min-k', required=False, default=8, type=int
)

args = parser.parse_args()


min_k = 8


with open(args.pickle, 'rb') as pkl:
    mate_pairs = pickle.load(pkl)

pairs_list = []
for _, r1, r2 in mate_pairs:
    pair = Pair(f=r1, r=r2) if r1.name.endswith('f') else Pair(f=r2, r=r1)
    if pair.f.comp == 'C':
        pair.f._comp_seq()
    if pair.r.comp != 'C':
        pair.r._comp_seq()
    pairs_list.append(pair)

output_fas = f"{args.output_dir}/concatenated_mate_pairs.fas"
with open(output_fas, 'w') as fas:
    output = ""
    for pair in pairs_list:
        output += f">{pair.f.name[:-1]}\n{pair.f.seq}{pair.r.seq}"
        print(output, file=fas)
        output = ""

einout = f"{args.output_dir}/concatenated_mate_pairs.einverted"
einfas = f"{args.output_dir}/concatenated_mate_pairs_tirs.fas"

subprocess.check_call(
    [
        'einverted', '-sequence', output_fas,
        '-outfile', einout, '-outseq', einfas,
        '-gap', '12', '-threshold', '50',
        '-match', '3', '-mismatch', '-4'
    ]
)
