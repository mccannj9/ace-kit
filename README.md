# ace-kit

Parsing and analysis of ace files using python

## Installation

Use the yaml file to create a conda environment

``` bash
conda env create -f ace-kit.yml
```

## Usage
To generate sequence logos and boundary plots for contigs in a cluster
using the default settings simply use:

``` python
python start.py -a [acefile] -o [output_directory]
```

This will give you the logos, plots, a fasta file with the sequence at
the inferred boundary for each contig and a pickled list of read pairs
where at least one of the mates is on a boundary _and_ the pair is not
split into separate clusters. This file, which ends with the extension
`.pkl` will be used as input for the next script.

``` python
python pairs_analysis [pickle_file]
```

The output of this script is still a bit of a work in progress, but will
eventually contain information about the presence of inverted repeats in
each read pair extracted from the previous analysis.
