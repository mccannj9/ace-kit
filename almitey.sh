#! /usr/bin/env bash

# $1 - the very top level cluster directory which contains seqclust
# $2 - the suffices of your paired reads as a string of length 2, such as "fr" or "12"

eval "$(conda shell.bash hook)"
conda activate ace-kit

clustering_dir="${1}/seqclust/clustering/clusters"
echo $clustering_dir

for d in ${clustering_dir}/dir_CL*
do
    echo "Working on cluster: ${d}"
    python almitey.py -i ${d} -o ${d}/almitey -s $2
done