# ace-kit

Parsing and analysis of ace files using python

## Installation

Use the yaml file to create a conda environment

``` bash
conda env create -f ace-kit.yml
```

## Usage

To run almitey on all of the clusters in a repeat explorer output simply
run the following command with the ```clustering_directory``` as the top
level directory of the repeat explorer output, i.e. the one containing the
seqclust subdirectory.

``` bash
python test_run.py -i [clustering_directory]
```

This will run the almitey analysis on every cluster in the output from the
RepeatExplorer and output a HTML file (almitey_report.html) in the same
folder where you find the your seqclust subdirectory. This report contains
information about almitey results for each cluster and allows navigation to
individual HTML reports, sequence logos and alignments for each cluster where
at least one repeat boundary was found

Use the following to see information about arguments to the script:

``` bash
python test_run.py --help
```

From this help menu you will see how you can change the default arguments when
running the analysis on your data.

If you want to run the analysis on an individual cluster you should either use
the ipython shell or write your own python script to import the Almitey class,
for example:

``` python
from almitey import Almitey

input_dir = "/home/user/clsutering/seqclust/clustering/clusters/dir_CL0045"
output_dir = "/home/user/your_choice"

runner = Almitey(input_dir, output_dir) # add other args to change defaults
results = runner.run() # add overwrite=True if necessary
```

By default this will not overwrite the output_dir directory, so you either have
to delete each time or use ```runner.run(overwrite=True)```
