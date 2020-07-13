#! /usr/bin/env python3

import os
import glob
import sys

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
        self.all = False

        if not(logfile):
            self.log_filename = f"{self.output_dir}/logfile.txt"

    def run(self, cluster):

        cluster_output_dict = {
            'cluster': '',
            'num_contigs': 0,
            'num_boundaries': 0,
            'avg_contig_length': 0,
            'avg_boundary_score': 0,
            'minor_path': ""
        }

        clname = os.path.basename(self.input_dir).split("_")[-1]
        cluster_output_dict['cluster'] = clname

        try:
            self.ace_filename = glob.glob(f"{self.input_dir}/*.ace")[0]

        except IndexError:
            print(f"No ace file found in {self.input_dir}", file=sys.stderr)
            return cluster_output_dict

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
            if len(sorted_contigs):
                cluster_output_dict['avg_contig_length'] = round(sum([
                    x.length for x in sorted_contigs
                ]) / len(sorted_contigs))

            nboundaries = sum([c.nboundaries for c in sorted_contigs])
            print(f"Total boundaries found: {nboundaries}", file=log)
            cluster_output_dict['num_boundaries'] = nboundaries

            if nboundaries:
                boundaries = []
                for c in sorted_contigs:
                    boundaries += c.boundaries
                boundaries.sort(key=lambda x: x.rate, reverse=True)
                cluster_output_dict['avg_boundary_score'] = round(sum([
                    x.rate for x in boundaries
                ]) / len(boundaries))

                for b in boundaries:
                    b.logo_path = os.path.basename(b.logo_path)
                dirname = os.path.basename(self.input_dir)

                if self.all:
                    cluster_output_dict['minor_path'] = f"{self.relative_loc}/{dirname}/almitey/almitey.html"

                with open(f"{self.output_dir}/almitey.html", 'w') as html:
                    clname = boundaries[0].contig.name.split("Contig")[0]
                    html_text = build_html_output(clname, boundaries)
                    print(html_text, file=html)

        return cluster_output_dict

    def run_on_all_clusters(self):
        self.all = True
        clusters = glob.glob(
            f"{self.input_dir}/{self.relative_loc}/dir_CL*"
        )
        clusters.sort(
            key=lambda c: int(os.path.basename(c).split("_")[-1][2:])
        )

        self.major_html_fn = f"{self.input_dir}/almitey_report.html"
        table_rows = []

        for cluster in clusters:
            print(cluster)
            self.input_dir = cluster
            self.output_dir = f"{cluster}/almitey"
            self.log_filename = f"{self.output_dir}/almitey_log.txt"
            result = self.run(cluster)

            if result['minor_path']:
                table_rows.append(major_row_template.safe_substitute(result))
            else:
                table_rows.append(major_row_none_template.safe_substitute(result))

        table_rows = "\n".join(table_rows)

        with open(self.major_html_fn, 'w') as major:
            html = major_html_template.safe_substitute({'table_rows': table_rows})
            print(html, file=major)
