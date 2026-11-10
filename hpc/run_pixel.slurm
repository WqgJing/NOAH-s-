#!/bin/bash
#SBATCH --job-name=noah_pixel
#SBATCH --output=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah/logs/pixel_%x_%j.out
#SBATCH --account=gts-rbras6-coda20 --partition=cpu-small -q inferno
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G --time=02:00:00
# Forward run of the PATCHED Noah-MP at one test pixel with our inputs and settings (NAME=BV|VR|KG, Y0, Y1):
#   ERA5 hourly forcing, Yuan LAI (DYNAMIC_VEG_OPTION 8), SMAP L4 soil moisture prescribed (PRESCRIBE_SOIL_MOISTURE 1),
#   Laio beta (BTR_OPTION 4, s* / s_w from the bank), L4 Campbell soil (SOIL_DATA_OPTION 4 + soil_spatial.nc).
#   Everything else = NCAR's shipped single-point namelist. FORCING 3600 s, model / output 1800 s, daily output files.
set -e
NAME=${NAME:?}; Y0=${Y0:-2018}; Y1=${Y1:-2019}
N=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah
PATCHED=$N/build/hrldas/hrldas
R=$N/runs/pixel_${NAME}_${Y0}_${Y1}
source /etc/profile.d/z00_lmod_pace.sh   # lmod is not initialised in non-login batch shells
source ~/.bashrc
module load gcc/12.3.0 mvapich2/2.3.7-1 netcdf-c/4.9.2-mva2-hdf5-1.14 netcdf-fortran/4.6.1-mva2-hdf5-1.14 || true
conda activate gw_viz

rm -rf $R; mkdir -p $R/forcing $R/out
SECONDS=0; python -u $N/scripts/make_pixel_forcing.py $NAME $Y0 $Y1 $R/forcing; echo "forcing: ${SECONDS}s"
source <(grep -E "^(LAIO_SSTAR|LAIO_SWILT|KDAY)=" $R/forcing/pixel_info.txt)

cd $R
cp $PATCHED/run/*.TBL .; ln -sf $PATCHED/run/hrldas.exe .
sed -e "s#^ *HRLDAS_SETUP_FILE *=.*# HRLDAS_SETUP_FILE = \"$R/forcing/hrldas_setup.nc\"#" \
    -e "s#^ *INDIR *=.*# INDIR = \"$R/forcing/\"#" -e "s#^ *OUTDIR *=.*# OUTDIR = \"$R/out/\"#" \
    -e "s#^ *START_YEAR *=.*# START_YEAR  = $Y0#" -e "s#^ *START_MONTH *=.*# START_MONTH = 01#" \
    -e "s#^ *START_DAY *=.*# START_DAY   = 01#" -e "s#^ *START_HOUR *=.*# START_HOUR  = 00#" -e "s#^ *KDAY *=.*# KDAY = $KDAY#" \
    -e "s#^ *FORCING_TIMESTEP *=.*# FORCING_TIMESTEP = 3600#" -e "s#^ *NOAH_TIMESTEP *=.*# NOAH_TIMESTEP    = 1800#" \
    -e "s#^ *OUTPUT_TIMESTEP *=.*# OUTPUT_TIMESTEP  = 1800#" -e "s#^ *SPLIT_OUTPUT_COUNT *=.*# SPLIT_OUTPUT_COUNT = 48#" \
    -e "s#^ *RESTART_FREQUENCY_HOURS *=.*# RESTART_FREQUENCY_HOURS = 8760#" \
    -e "s#^ *DYNAMIC_VEG_OPTION *=.*# DYNAMIC_VEG_OPTION                = 8#" \
    -e "s#^ *BTR_OPTION *=.*# BTR_OPTION                        = 4#" \
    -e "s#^ *SOIL_DATA_OPTION *=.*# SOIL_DATA_OPTION                  = 4#" \
    -e "s#^ *CROP_OPTION *=.*# CROP_OPTION                       = 0#" \
    $N/src/hrldas/hrldas/run/examples/single_point/namelist.hrldas.single_point > namelist.tmp
sed -i "0,/^\/ *$/s##  SPATIAL_FILENAME = \"$R/forcing/soil_spatial.nc\"\n  PRESCRIBE_SOIL_MOISTURE = 1\n  LAIO_SSTAR = $LAIO_SSTAR\n  LAIO_SWILT = $LAIO_SWILT\n/#" namelist.tmp
mv namelist.tmp namelist.hrldas
echo "== namelist vs NCAR example:"; diff $N/src/hrldas/hrldas/run/examples/single_point/namelist.hrldas.single_point namelist.hrldas | grep "^>" || true
SECONDS=0; ./hrldas.exe > hrldas.log 2>&1 || { echo "hrldas.exe FAILED"; tail -30 hrldas.log; exit 1; }
echo "hrldas.exe: ${SECONDS}s; $(grep 'SFLXtot' hrldas.log | tail -1 | sed 's/  */ /g')"
grep -i -E "error|warning|stop" hrldas.log | sort | uniq -c | head -10 || true
echo "output files: $(ls out | grep -c LDASOUT)"
rm -f $R/forcing/*.LDASIN_DOMAIN1
echo RUN_OK
