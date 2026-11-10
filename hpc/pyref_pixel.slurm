#!/bin/bash
#SBATCH --job-name=pyref
#SBATCH --output=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah/logs/pyref_%x_%j.out
#SBATCH --account=gts-rbras6-coda20 --partition=cpu-small -q inferno
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G --time=01:00:00
# REFERENCE for the Python port (noah/pyport): patched Fortran Noah-MP, our full configuration (as hpc/run_pixel.slurm),
# at pixel NAME for a snow-free window: cold start from the setup file at START (default 2018-06-01), KDAY days (default 122).
# Keeps forcing/forcing_arrays.npz (inputs for Python) and out/*.LDASOUT_DOMAIN1 (reference outputs); LDASIN files are deleted.
set -e
NAME=${NAME:?}; START=${START:-2018-06-01}; KDAY=${KDAY:-122}
N=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah
PATCHED=$N/build/hrldas/hrldas
R=$N/runs/pyref_${NAME}
Y=${START:0:4}
source /etc/profile.d/z00_lmod_pace.sh; source ~/.bashrc
module load gcc/12.3.0 mvapich2/2.3.7-1 netcdf-c/4.9.2-mva2-hdf5-1.14 netcdf-fortran/4.6.1-mva2-hdf5-1.14 || true
conda activate gw_viz
rm -rf $R; mkdir -p $R/forcing $R/out
python -u $N/scripts/make_pixel_forcing.py $NAME $Y $Y $R/forcing | grep -E "setup check|wrote|L4 root"
source <(grep -E "^(LAIO_SSTAR|LAIO_SWILT)=" $R/forcing/pixel_info.txt)
cd $R; cp $PATCHED/run/*.TBL .; ln -sf $PATCHED/run/hrldas.exe .
sed -e "s#^ *HRLDAS_SETUP_FILE *=.*# HRLDAS_SETUP_FILE = \"$R/forcing/hrldas_setup.nc\"#" \
    -e "s#^ *INDIR *=.*# INDIR = \"$R/forcing/\"#" -e "s#^ *OUTDIR *=.*# OUTDIR = \"$R/out/\"#" \
    -e "s#^ *START_YEAR *=.*# START_YEAR  = $Y#" -e "s#^ *START_MONTH *=.*# START_MONTH = ${START:5:2}#" \
    -e "s#^ *START_DAY *=.*# START_DAY   = ${START:8:2}#" -e "s#^ *START_HOUR *=.*# START_HOUR  = 00#" -e "s#^ *KDAY *=.*# KDAY = $KDAY#" \
    -e "s#^ *FORCING_TIMESTEP *=.*# FORCING_TIMESTEP = 3600#" -e "s#^ *NOAH_TIMESTEP *=.*# NOAH_TIMESTEP    = 1800#" \
    -e "s#^ *OUTPUT_TIMESTEP *=.*# OUTPUT_TIMESTEP  = 1800#" -e "s#^ *SPLIT_OUTPUT_COUNT *=.*# SPLIT_OUTPUT_COUNT = 48#" \
    -e "s#^ *RESTART_FREQUENCY_HOURS *=.*# RESTART_FREQUENCY_HOURS = 8760#" \
    -e "s#^ *DYNAMIC_VEG_OPTION *=.*# DYNAMIC_VEG_OPTION                = 8#" -e "s#^ *BTR_OPTION *=.*# BTR_OPTION                        = 4#" \
    -e "s#^ *SOIL_DATA_OPTION *=.*# SOIL_DATA_OPTION                  = 4#" -e "s#^ *CROP_OPTION *=.*# CROP_OPTION                       = 0#" \
    $N/src/hrldas/hrldas/run/examples/single_point/namelist.hrldas.single_point > namelist.hrldas
sed -i "0,/^\/ *$/s##  SPATIAL_FILENAME = \"$R/forcing/soil_spatial.nc\"\n  PRESCRIBE_SOIL_MOISTURE = 1\n  LAIO_SSTAR = $LAIO_SSTAR\n  LAIO_SWILT = $LAIO_SWILT\n/#" namelist.hrldas
./hrldas.exe > hrldas.log 2>&1 || { echo "hrldas.exe FAILED"; head -5 hrldas.log; exit 1; }
NOUT=$(ls out | grep -c LDASOUT || true); echo "output files: $NOUT (expected $KDAY)"
rm -f $R/forcing/*.LDASIN_DOMAIN1
[ "$NOUT" -ge "$KDAY" ] && echo PYREF_OK || { echo PYREF_FAILED; exit 1; }
