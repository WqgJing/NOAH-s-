#!/bin/bash
#SBATCH --job-name=pyrefv
#SBATCH --output=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah/logs/pyrefv_%x_%j.out
#SBATCH --account=gts-rbras6-coda20 --partition=cpu-small -q inferno
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:30:00
# Variant of a Python-port reference run (hpc/pyref_pixel.slurm): same inputs as runs/pyref_${NAME}, other Laio s*, s_w
# and optional nuisance scales, written to runs/pyref_${NAME}${TAG}. Used to test the BTRAN < 1 branch of noah/pyport.
#   sbatch --export=ALL,NAME=VR,TAG=s,SSTAR=0.25,SWILT=0.10 hpc/pyref_variant.slurm
#   sbatch --export=ALL,NAME=KG,TAG=sn,SSTAR=0.25,SWILT=0.10,SCALES="1.5 0.7 1.3" hpc/pyref_variant.slurm
set -e
NAME=${NAME:?}; TAG=${TAG:?}; SSTAR=${SSTAR:?}; SWILT=${SWILT:?}; SCALES=${SCALES:-}
N=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah
SRC=$N/runs/pyref_${NAME}; R=$N/runs/pyref_${NAME}${TAG}
source /etc/profile.d/z00_lmod_pace.sh; source ~/.bashrc
module load gcc/12.3.0 mvapich2/2.3.7-1 netcdf-c/4.9.2-mva2-hdf5-1.14 netcdf-fortran/4.6.1-mva2-hdf5-1.14 || true
conda activate gw_viz
rm -rf $R; mkdir -p $R
python -u $N/scripts/make_pyref_variant.py $SRC $R $SSTAR $SWILT $SCALES
cd $R; ln -sf $N/build/hrldas/hrldas/run/hrldas.exe .
grep -E "LAIO|NUIS|SETUP|SPATIAL|INDIR|OUTDIR" namelist.hrldas
./hrldas.exe > hrldas.log 2>&1 || { echo "hrldas.exe FAILED"; head -5 hrldas.log; exit 1; }
NOUT=$(ls out | grep -c LDASOUT || true); NREF=$(ls $SRC/out | grep -c LDASOUT || true)
echo "output files: $NOUT (reference run $NREF)"
rm -f $R/forcing/*.LDASIN_DOMAIN1
[ "$NOUT" -ge "$NREF" ] && echo PYREF_OK || { echo PYREF_FAILED; tail -5 hrldas.log; exit 1; }
