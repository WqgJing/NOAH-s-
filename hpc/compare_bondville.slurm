#!/bin/bash
#SBATCH --job-name=bondville_cmp
#SBATCH --output=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah/logs/bondville_cmp_%j.out
#SBATCH --account=gts-rbras6-coda20 --partition=cpu-small -q inferno
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:20:00
# Compare our Bondville runs (shipped + tutorial options) with the NCAR tutorial's run log and LH screenshot.
set -e
cd /storage/project/r-rbras6-0/jjason3/proposal_drought/noah; source ~/.bashrc; conda activate gw_viz
python -u scripts/compare_bondville.py
