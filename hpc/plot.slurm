#!/bin/bash
#SBATCH --job-name=noah_plot
#SBATCH --output=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah/logs/plot_%j.out
#SBATCH --account=gts-rbras6-coda20 --partition=cpu-small -q inferno
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:20:00
# Run a plotting script from noah/scripts: sbatch --export=ALL,SCRIPT=plot_bondville_trad.py hpc/plot.slurm
set -e
cd /storage/project/r-rbras6-0/jjason3/proposal_drought/noah; source ~/.bashrc; conda activate gw_viz
python -u scripts/${SCRIPT:?set SCRIPT}
