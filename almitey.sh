#! /usr/bin/env bash

eval "$(conda shell.bash hook)"
conda activate ace-kit

clustering_dir="${1}/seqclust/clustering/clusters"
echo $clustering_dir

for d in ${clustering_dir}/dir_CL*
do
    echo "Working on cluster: ${d}"
    echo "python almitey.py -i ${d} -o ${d}/almitey"
done