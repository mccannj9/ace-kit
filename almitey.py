#! /usr/bin/env python3

import os
import glob
import sys
import argparse

from kit.finder import SwitchpointFinder
from kit.blast import set_result_orientation, quick_blastn, parse_blast_output
from kit.utils import find_all_boundary_reads, pair_boundary_reads, pairs_with_correct_orient
from kit.html import build_html_output


class Almitey(object):
    def __init__(
        self, input_dir, output_dir, window_size=7, min_depth=10,
        min_read_prop=0.01, logfile=None, suffices="fr"
    ):

        self.input_dir = input_dir
        try:
            self.ace_filename = glob.glob(f"{args.input_dir}/*.ace")[0]

        except IndexError:
            print(f"No ace file found in {self.input_dir}", file=sys.stderr)
            sys.exit(0)

        self.output_dir = output_dir
        self.window_size = window_size
        self.min_depth = min_depth
        self.min_read_prop = min_read_prop
        self.log_filename = logfile
        self.suffices = suffices

    def run(self, cluster):
        try:
            os.mkdir(self.output_dir)

        except FileExistsError:
            print(f"{self.output_dir} exists already. Continuing...")

        with open(self.log_filename, 'w') as log:
            finder = SwitchpointFinder(
                self.ace_filename, outdir=self.output_dir, window_size=self.window_size,
                min_depth=self.min_depth, min_read_prop=self.min_read_prop
            )
            contigs, all_reads = finder.fit()
            sorted_contigs = sorted(
                contigs, key=lambda c: (c.nboundaries, c.boundary_rate_sum), reverse=True
            )
            # remove any contigs with no inferred boundaries
            sorted_contigs[:] = [x for x in sorted_contigs if len(x.boundaries)]

            nboundaries = sum([c.nboundaries for c in sorted_contigs])
            print(f"Total boundaries found: {nboundaries}", file=log)

            if not(nboundaries):
                return

            else:
                boundaries = []
                for c in sorted_contigs:
                    boundaries += c.boundaries
                boundaries.sort(key=lambda x: x.rate, reverse=True)

            with open(f"{self.output_dir}/almitey.html", 'w') as html:
                clname = boundaries[0].contig.name.split("Contig")[0]
                html_text = build_html_output(clname, boundaries)
                print(html_text, file=html)


    def run_on_all_clusters(self):
        clusters = glob.glob(
            f"{self.input_dir}/seqclust/clustering/clusters/dir_CL*"
        )

        for cluster in clusters:
            self.output_dir = f"{cluster}/almitey"
            self.log_filename = f"{self.output_dir}/almitey_log.txt"
            self.run(cluster)


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

try:
    acefile = glob.glob(f"{args.input_dir}/*.ace")[0]

except IndexError:
    print("No ace file found", file=sys.stderr)
    sys.exit(0)

try:
    os.mkdir(args.output_dir)

except FileExistsError:
    pass

if args.log_file:
    log = open(f"{args.output_dir}/{args.log_file}", 'w')

else:
    log = sys.stderr

suffix_1 = args.read_suffices[0]
suffix_2 = args.read_suffices[-1]

keyword_args = {
    'window_size': args.window_size,
    'min_depth': args.min_depth,
    'min_read_prop': args.min_read_prop
}

finder = SwitchpointFinder(acefile, args.output_dir, **keyword_args)
contigs, all_reads = finder.fit()

sorted_contigs = sorted(
    contigs, key=lambda c: (c.nboundaries, c.boundary_rate_sum), reverse=True
)
# remove any contigs with no inferred boundaries
sorted_contigs[:] = [x for x in sorted_contigs if len(x.boundaries)]

nboundaries_found = sum([c.nboundaries for c in sorted_contigs])

if not(nboundaries_found):
    print("Not enough boundaries found in contigs")
    # run was still a success, just no boundaries ;)
    print("No boundaries found", file=log)
    sys.exit(0)

# at_least_two = True if nboundaries_found > 0 else False
# check if two boundaries were found in one contig
two = True if sorted_contigs[0].nboundaries == 2 else False

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
    orient_boundaries = True

else:
    print("Method for two boundaries in different contigs not ready yet", file=log)
    orient_boundaries = False

boundaries = []
for c in sorted_contigs:
    boundaries += c.boundaries

if orient_boundaries:
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

boundary_reads = find_all_boundary_reads(boundaries)
paired_boundary_reads = pair_boundary_reads(
    boundary_reads, all_reads, suffix_1=suffix_1, suffix_2=suffix_2
)
oriented_pairs = pairs_with_correct_orient(paired_boundary_reads)

with open(f"{finder.outdir}/potential_boundary_pairs.fas", 'w') as fasta:
    for k in oriented_pairs:
        l, r = oriented_pairs[k]
        lout = f'>{l.name} {l.ctg} {l.comp} {int(l.side)}\n{l.seq.replace("*", "")}'
        rout = f'>{r.name} {r.ctg} {r.comp} {int(r.side)}\n{r.seq.replace("*", "")}'
        print(lout, file=fasta)
        print(rout, file=fasta)

try:
    perc_oriented = len(oriented_pairs) / len(paired_boundary_reads)

except ZeroDivisionError:
    print("No boundary pairs found.")

print(f"Oriented Pairs: {len(oriented_pairs)}", file=log)
print(f"Total Paired with 1 Boundary: {len(paired_boundary_reads)}", file=log)

with open(f"{finder.outdir}/almitey.html", 'w') as html:
    clname = boundaries[0].contig.name.split("Contig")[0]
    html_text = build_html_output(clname, boundaries)
    print(html_text, file=html)
