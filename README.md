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

Use the following to see information about arguments to the script:

```python
python start.py --help
```

``` python
python pairs_analysis -p [pickle_file] -o [output_dir]
```

The output of this script is currently the results of running the program
einverted from the emboss package. This includes an alignment file, similar
to what one gets when running a blast with the standard blast output format.
There is also a fasta file with the extracted regions of the reads which
have inverted repeats. Arguments for parameters that can be tweaked for the
einverted program will be added, as of now they are fixed to defaults.
