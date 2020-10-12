
from string import Template

minor_table_row = """
                <div class='table-body-row'>
                    <div class="table-body-cell">
                        $contig_name
                    </div>
                    <div class="table-body-cell">
                        $pos
                    </div>
                    <div class="table-body-cell">
                        $depth
                    </div>
                    <div class="table-body-cell">
                        $rate
                    </div>
                    <div class="table-body-cell">
                        $side
                    </div>
                    <div class="table-logo-cell">
                        <a href="$logo_path"><img src="$logo_path" alt="Logo Plot" width=400></a>
                    </div>
                    <div class="table-seq-cell">
                        $seq
                    </div>
                </div>
"""

minor_html = """
<html>
    <head>
        <title>Almitey Report $cluster_name</title>
        <link rel="stylesheet" type="text/css" href="results.css">
    </head>
    <style>
        #table{width:100%;display:table;border-collapse:collapse}#table-header{display:table-header-group;font-weight:700;font-size:25px}#table-body{display:table-row-group}.table-header-cell{display:table-cell;padding:10px;text-align:center;border-bottom:1px solid #000}.table-body-row{display:table-row;border-bottom:1px solid grey}.table-body-cell{display:table-cell;text-align:center;vertical-align:middle}.table-logo-cell{display:table-cell;vertical-align:middle;text-align:center}.table-seq-cell{display:table-cell;font-family:monospace;text-align:center;vertical-align:middle}
    </style>
    <body>
        <h1 class="cluster-result-header">Almitey Report $cluster_name</h1>
        </table>
        <div id="table">
            <div id="table-header">
                <div class="table-header-cell">Contig</div>
                <div class="table-header-cell">Position</div>
                <div class="table-header-cell">Depth</div>
                <div class="table-header-cell">Score</div>
                <div class="table-header-cell">Contig Side</div>
                <div class="table-header-cell">Logo</div>
                <div class="table-header-cell">Sequence</div>
            </div>
            <div id="table-body">
                $table_rows
            </div>
        </div>
    </body>
</html>
"""

major_table_row = """
                <div class='table-body-row'>
                    <div class="table-body-cell">
                        <a href="$minor_path">$cluster</a>
                    </div>
                    <div class="table-body-cell">
                        $num_contigs
                    </div>
                    <div class="table-body-cell">
                        $num_boundaries
                    </div>
                    <div class="table-body-cell">
                        $avg_contig_length
                    </div>
                    <div class="table-body-cell">
                        $avg_boundary_score
                    </div>
                    <div class="table-body-cell">
                        $alignment_path
                        <!--<a href="$alignment_path">alignment</a>-->
                    </div>
                    <div class="table-body-cell">
                        $TIR
                    </div>
                </div>
"""

major_table_row_none = """
                <div class='table-body-row'>
                    <div class="table-body-cell">
                        $cluster
                    </div>
                    <div class="table-body-cell">
                        $num_contigs
                    </div>
                    <div class="table-body-cell">
                        $num_boundaries
                    </div>
                    <div class="table-body-cell">
                        $avg_contig_length
                    </div>
                    <div class="table-body-cell">
                        $avg_boundary_score
                    </div>
                    <div class="table-body-cell">
                        alignment
                    </div>
                    <div class="table-body-cell">
                        NO
                    </div>
                </div>
"""

major_html = """
<html>
    <head>
        <title>Complete by Cluster Almitey Report</title>
        <link rel="stylesheet" type="text/css" href="results.css">
    </head>
    <style>
        #table{width:100%;display:table;border-collapse:collapse}#table-header{display:table-header-group;font-weight:700;font-size:25px}#table-body{display:table-row-group}.table-header-cell{display:table-cell;padding:10px;text-align:center;border-bottom:1px solid #000}.table-body-row{display:table-row;border-bottom:1px solid grey}.table-body-cell{display:table-cell;text-align:center;vertical-align:middle}.table-logo-cell{display:table-cell;vertical-align:middle;text-align:center}.table-seq-cell{display:table-cell;font-family:monospace;text-align:center;vertical-align:middle}
    </style>
    <body>
        <h1 class="cluster-result-header">Complete by Cluster Almitey Report</h1>
        </table>
        <div id="table">
            <div id="table-header">
                <div class="table-header-cell">Cluster</div>
                <div class="table-header-cell">No. Contigs</div>
                <div class="table-header-cell">No. Boundaries</div>
                <div class="table-header-cell">Avg. Contig Length</div>
                <div class="table-header-cell">Avg. Boundary Score</div>
                <div class="table-header-cell">Alignment</div>
                <div class="table-header-cell">TIR?</div>
            </div>
            <div id="table-body">
                $table_rows
            </div>
        </div>
    </body>
</html>
"""

minor_row_template = Template(minor_table_row)
minor_html_template = Template(minor_html)
major_row_template = Template(major_table_row)
major_row_none_template = Template(major_table_row_none)
major_html_template = Template(major_html)

def build_html_output(cluster_name: str, boundaries: list) -> str:
    table_rows = []
    for b in boundaries:
        table_rows.append(b.table_row_template().strip('\n'))
    table_rows = "\n".join(table_rows)
    html = minor_html_template.safe_substitute({
        'cluster_name': cluster_name,
        'table_rows': table_rows
    })

    return html
