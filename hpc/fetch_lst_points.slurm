#!/bin/bash
#SBATCH --job-name=lstpoints
#SBATCH --output=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah/logs/lstpoints_%j.out
#SBATCH --account=gts-rbras6-coda20 --partition=cpu-small -q inferno
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=02:00:00
# MODIS Aqua/Terra daily LST at the Noah-MP test pixels (GEE, same service account as the L4 extractors).
set -e
cd /storage/project/r-rbras6-0/jjason3/proposal_drought/noah; mkdir -p data; source ~/.bashrc; conda activate gw_viz
python -u scripts/fetch_lst_points.py
