
import os
import subprocess
from dataclasses import dataclass, fields

from kit.utils import sign

blast_path = os.environ['blast_path'] if 'blast_path' in os.environ else 'blastn'
output_dir = os.environ['output_dir'] if 'output_dir' in os.environ else './'


@dataclass
class BlastResult:
    query: str
    subject: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    orientation: int = 0

    def __eq__(self, other):
        q = self.query == other.query
        s = self.subject == other.subject
        return q and s



def quick_blastn(
    query:str=None, subject:str=None, out:str=None,
    outfmt:str='6', add_args=['-task', 'blastn-short']
):
    if not(query):
        print("please provide at least a query file")
        return 1

    if not(subject):
        subject = query

    if not(out):
        out = f"{output_dir}/boundary_blast_output.txt"

    try:
        args = [
            blast_path, '-query', query, '-subject',
            subject, '-out', out, '-outfmt', outfmt,
        ] + add_args
        subprocess.check_call(args)

    except subprocess.CalledProcessError as e:
        call = " ".join(args)
        print(f"The following call returned the error code = {e}:\n{call}")
        return e

    return out


def parse_blast_output(input_file):
    blast_results = []
    with open(input_file) as result_file:
        for line in result_file:
            line = line.strip().split()
            res = BlastResult(*[
                fields(BlastResult)[x].type(y) for x, y in enumerate(line)
            ])

            if res in blast_results:
                other = blast_results[blast_results.index(res)]
                if res.evalue < other.evalue:
                    blast_results[blast_results.index(res)] = res
                continue

            blast_results.append(res)

    return blast_results


def set_result_orientation(br: BlastResult):
    qdiff = br.qstart - br.qend
    sdiff = br.sstart - br.send

    if sign(qdiff) == sign(sdiff):
        br.orientation = 1
    else:
        br.orientation = -1

    return br.orientation


def get_blast_hits_with_orientation(fasta:str, out:str=f"{output_dir}/boundary_blast_output.txt"):
    retval = quick_blastn(fasta, out=out)
    blast_res = parse_blast_output(retval)

    for hit in blast_res:
        set_blast_result_orientation(hit)

    return blast_res
