#! /usr/bin/env python3

import os
import glob
import sys
import argparse

from kit.finder import SwitchpointFinder
from kit.blast import set_result_orientation, quick_blastn, parse_blast_output
from kit.utils import find_all_boundary_reads, pair_boundary_reads, pairs_with_correct_orient

from kit.html import build_html_output, major_html_template
from kit.html import major_row_template, major_row_none_template

class Almitey(object):
    def __init__(
        self, input_dir, output_dir, window_size=7, min_depth=10,
        min_read_prop=0.01, logfile=None, suffices="fr"
    ):

        self.input_dir = input_dir
        self.output_dir = output_dir
        self.window_size = window_size
        self.min_depth = min_depth
        self.min_read_prop = min_read_prop
        self.log_filename = logfile
        self.suffices = suffices
        self.relative_loc = "seqclust/clustering/clusters"

    def run(self, cluster, all=False):

        cluster_output_dict = {
            'cluster': '',
            'num_contigs': 0,
            'num_boundaries': 0,
            'avg_contig_len': 0,
            'avg_boundary_score': 0,
            'minor_path': ""
        }

        clname = os.path.basename(cluster).split("_")[-1]
        cluster_output_dict['cluster'] = clname

        try:
            self.ace_filename = glob.glob(f"{cluster}/*.ace")[0]

        except IndexError:
            print(f"No ace file found in {self.input_dir}", file=sys.stderr)
            return major_row_none_template.safe_substitute(cluster_output_dict)

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
            cluster_output_dict['num_contigs'] = len(contigs)
            sorted_contigs = sorted(
                contigs, key=lambda c: (c.nboundaries, c.boundary_rate_sum), reverse=True
            )
            # remove any contigs with no inferred boundaries
            sorted_contigs[:] = [x for x in sorted_contigs if len(x.boundaries)]

            nboundaries = sum([c.nboundaries for c in sorted_contigs])
            print(f"Total boundaries found: {nboundaries}", file=log)
            cluster_output_dict['num_boundaries'] = nboundaries

            if nboundaries:
                boundaries = []
                for c in sorted_contigs:
                    boundaries += c.boundaries
                boundaries.sort(key=lambda x: x.rate, reverse=True)

                if all:
                    # setting up paths
                    for b in boundaries:
                        b.logo_path = os.path.basename(b.logo_path)
                    dirname = os.path.basename(cluster)
                    cluster_output_dict['minor_path'] = f"{self.relative_loc}/{dirname}/almitey/almitey.html"
                
                with open(f"{self.output_dir}/almitey.html", 'w') as html:
                    clname = boundaries[0].contig.name.split("Contig")[0]
                    html_text = build_html_output(clname, boundaries)
                    print(html_text, file=html)

        return major_row_template.safe_substitute(cluster_output_dict)


    def run_on_all_clusters(self):
        clusters = glob.glob(
            f"{self.input_dir}/{self.relative_loc}/dir_CL*"
        )
        clusters.sort(
            key=lambda c: int(os.path.basename(c).split("_")[-1][2:])
        )

        self.major_html_fn = f"{self.input_dir}/almitey.html"
        table_rows = []

        for cluster in clusters:
            print(cluster)
            self.output_dir = f"{cluster}/almitey"
            self.log_filename = f"{self.output_dir}/almitey_log.txt"
            html = self.run(cluster, all=True)
            table_rows.append(html)

        table_rows = "\n".join(table_rows)

        with open(self.major_html_fn, 'w') as major:
            html = major_html_template.safe_substitute({'table_rows': table_rows})
            print(html, file=major)
