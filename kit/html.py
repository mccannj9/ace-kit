
cluster_name = None

contig_number = None
path_boundary_plot = None
path_left_logo = None
path_right_logo = None

table_row = f"""
                <div class='table-body-row'>
                    <div class="table-body-cell">
                        {contig_number}
                    </div>
                    <div class="table-body-cell">
                        <img src="{path_boundary_plot}" alt="Masking Plot" width=400>
                    </div>
                    <div class="table-logo-cell">
                        <img src="{path_left_logo}" alt="Masking Plot" width=400>
                        <img src="{path_right_logo}" alt="Masking Plot" width=400>
                    </div>
                </div>

"""

table_rows = None

html = f"""
<html>
    <head>
        <title>Almitey Report {cluster_name}</title>
        <link rel="stylesheet" type="text/css" href="results.css">
    </head>
    <body>
        <h1 class="cluster-result-header">Almitey Report {cluster_name}</h1>
        </table>
        <div id="table">
            <div id="table-header">
                <div class="table-header-cell">Contig Number</div>
                <div class="table-header-cell">Masking Plot</div>
                <div class="table-header-cell">Logos</div>
            </div>
            <div id="table-body">
                {table_rows}
            </div>
        </div>
    </body>
</html>

"""

