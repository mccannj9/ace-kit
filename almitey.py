#! /usr/bin/env python3

import os
import glob
import sys
import argparse
import pickle

from kit.finder import SwitchpointFinder
from kit.blast import set_result_orientation, quick_blastn, parse_blast_output

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

acefile = glob.glob(f"{args.input_dir}/*.ace")[0]

finder = SwitchpointFinder(acefile, args.output_dir, **keyword_args)
contigs = finder.fit_new()

sorted_contigs = sorted(
    contigs, key=lambda c: (c.nboundaries, c.boundary_rate_sum), reverse=True
)
# remove any contigs with no inferred boundaries
sorted_contigs[:] = [x for x in sorted_contigs if len(x.boundaries)]

nboundaries_found = sum([c.nboundaries for c in sorted_contigs])
at_least_two = True if nboundaries_found > 0 else False
# check if two boundaries were found in one contig
two = True if sorted_contigs[0].nboundaries == 2 else False

if not(at_least_two):
    print("Not enough boundaries found in contigs")
    # run was still a success, just no boundaries ;)
    sys.exit(0)

# building database from top contig with two boundaries
if two:
    reference = sorted_contigs[0]
    outname = f"{args.output_dir}/top_boundaries_db.fas"
    print(f"Two boundaries found in {reference.name}, using as database > {outname}")

    with open(outname, 'w') as fasta:
        for b in reference.boundaries:
            b.orient = b.side
            print(b.boundary_seq_as_fasta(), file=fasta)

    ref_l, ref_r = reference.boundaries

else:
    print(f"Method for two boundaries in different contigs not ready yet")
    sys.exit(0)

boundaries = []
for c in sorted_contigs:
    boundaries += c.boundaries

query = f"{finder.outdir}/boundaries_from_contigs.fas"
subject = outname
orient_out = f"{finder.outdir}/orientation_blast.txt"
quick_blastn(query, subject, orient_out)

blast_results = parse_blast_output(orient_out)
for res in blast_results:
    set_result_orientation(res)

# don't look at first two boundaries, they are the reference
boundary_blasts = {}
for b in boundaries[2:]:
    boundary_blasts[b.name] = sorted(
        [x for x in blast_results if b.name == x.query], key=lambda b: b.subject
    )
    orients = [x.orientation for x in boundary_blasts[b.name]]
    if orients == [ref_l.orient, ref_r.orient]:
        print('left side')
        b.orient = ref_l.orient
    else:
        print('right side')
        b.orient = ref_r.orient


with open(f"{finder.outdir}/boundaries_left.fas", 'w') as fasta:
    for b in boundaries:
        if b.orient == 1.0:
            print(b.boundary_seq_as_fasta(), file=fasta)

with open(f"{finder.outdir}/boundaries_right.fas", 'w') as fasta:
    for b in boundaries:
        if b.orient == -1.0:
            print(b.boundary_seq_as_fasta(), file=fasta)


# paired reads analysis

subject_prefix = f"{finder.outdir}/boundaries"
reads_fn = f"{args.input_dir}/reads.fas"

print("Blasting all reads against reference boundaries..")
quick_blastn(
    query=reads_fn, subject=outname, out=f"{finder.outdir}/all_reads_blast.txt"
)
print("Finished.")

blast_results = parse_blast_output(f"{finder.outdir}/all_reads_blast.txt")

id_dict = {}
for br in blast_results:
    n = br.query[:-1]
    if n not in id_dict:
        id_dict[n] = 1
    else:
        id_dict[n] += 1

potential_pairs = set([x for x in id_dict if id_dict[x] == 4])
