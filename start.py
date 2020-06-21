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
acefile = finder.acefile

results_dict, all_reads = finder.fit()
results_list = sorted(
    results_dict.values(),
    key=lambda element: (element.n, element.max_mag_derivative()),
    reverse=True
)

# checking if first contig in sorted results has 2 boundaries
two = True if results_list[0].n == 2 else False
at_least_one = True if results_list[0].n > 0 else False

if not(at_least_one):
    print("No boundaries found in contigs")
    # run was still a success, just no boundaries ;)
    sys.exit(0)

# drop contig results with no boundaries found
# possible because a contig shows up here if coverage is over cov threshold
results_list[:] = [x for x in results_list if x.n]

# deal with these two cases, two boundaries in one ctg is easiest
if two:
    db_result = results_list[0]
    outfilename = f"{args.output_dir}/top_boundaries_db.fas"
    print(f"Two boundaries found in one contig, using as database > {outfilename}")
    db_result.write_contig_boundaries_as_fasta(outfilename)


else:
    # get 2 contigs with boundaries on opposite sides
    pass

query = f"{finder.outdir}/boundaries_from_contigs.fas"
subject = outfilename
orient_out = f"{finder.outdir}/orientation_blast.txt"
quick_blastn(query, subject, orient_out)
blast_results = parse_blast_output(orient_out)
for res in blast_results:
    set_blast_result_orientation(res)
br_sorted = sorted(
    blast_results,
    key=lambda element: (element.subject, element.orientation)
)

keep_results = []
for res in results_list:
    contig_name = res.contig.name
    for br in br_sorted:
        if br.query.split("_")[0] == contig_name:
            res.blast_results.append(br)
    if len(res.blast_results) == res.n * 2:
        keep_results.append(res)
results = finder.results

boundary_read_prefixes = set([y.name[:-1] for x,y in all_reads.items() if y.boundary])

mate_pairs = []
for prefix in boundary_read_prefixes:
    f = f"{prefix}f"
    r = f"{prefix}r"
    if (f in all_reads) and (r in all_reads):
        rf = all_reads[f]
        rr = all_reads[r]
        remove_gaps_and_adjust_boundary(rf)
        remove_gaps_and_adjust_boundary(rr)
        mate_pairs.append(
            (rf.boundary * rr.boundary, rf, rr)
        )

boundary_mate_pairs = [(x,y) for w,x,y in mate_pairs if w]

with open(f'{args.output_dir}/mates_with_boundary_info.pkl', 'wb') as pkl:
    pickle.dump(mate_pairs, pkl)

fasta = f"{finder.outdir}/boundaries_from_contigs.fas"
blast_output = f"{finder.outdir}/boundary_blast_output.txt"
blast_hits = get_blast_hits_with_orientation(fasta, out=blast_output)
