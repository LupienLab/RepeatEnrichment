#!/bin/bash
#SBATCH -p long
#SBATCH --job-name pipeline
#SBATCH -t 7-00:00:00
#SBATCH -o log.out
#SBATCH -e log.err

snakemake -j 2 -c "sbatch {cluster.params}" -u config_slurm/slurm.yaml --latency-wait 60 --nolock
