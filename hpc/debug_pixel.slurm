#!/bin/bash
#SBATCH --job-name=noah_debug
#SBATCH --output=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah/logs/debug_%j.out
#SBATCH --account=gts-rbras6-coda20 --partition=cpu-small -q inferno
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G --time=01:00:00
# Debug a pixel run that stops after the first step: 3-day runs (2018-01-01..04) with the full configuration and with one switch
# turned back to the stock setting each time; prints the exit code and the end of each log. Forcing for 2018 only, kept.
NAME=${NAME:-KG}
N=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah
PATCHED=$N/build/hrldas/hrldas
R=$N/runs/debug_${NAME}
source /etc/profile.d/z00_lmod_pace.sh; source ~/.bashrc
module load gcc/12.3.0 mvapich2/2.3.7-1 netcdf-c/4.9.2-mva2-hdf5-1.14 netcdf-fortran/4.6.1-mva2-hdf5-1.14 || true
conda activate gw_viz
export GFORTRAN_ERROR_BACKTRACE=1
rm -rf $R; mkdir -p $R/forcing
python -u $N/scripts/make_pixel_forcing.py $NAME 2018 2018 $R/forcing | grep -E "setup check|wrote"
source <(grep -E "^(LAIO_SSTAR|LAIO_SWILT)=" $R/forcing/pixel_info.txt)
case_run () {   # $1 name, then sed expressions applied on top of the full configuration
  local C=$1; shift; local D=$R/$C; mkdir -p $D/out; cd $D; cp $PATCHED/run/*.TBL .; ln -sf $PATCHED/run/hrldas.exe .
  sed -e "s#^ *HRLDAS_SETUP_FILE *=.*# HRLDAS_SETUP_FILE = \"$R/forcing/hrldas_setup.nc\"#" \
      -e "s#^ *INDIR *=.*# INDIR = \"$R/forcing/\"#" -e "s#^ *OUTDIR *=.*# OUTDIR = \"$D/out/\"#" \
      -e "s#^ *START_YEAR *=.*# START_YEAR  = 2018#" -e "s#^ *START_DAY *=.*# START_DAY   = 01#" -e "s#^ *KDAY *=.*# KDAY = 3#" \
      -e "s#^ *FORCING_TIMESTEP *=.*# FORCING_TIMESTEP = 3600#" -e "s#^ *SPLIT_OUTPUT_COUNT *=.*# SPLIT_OUTPUT_COUNT = 48#" \
      -e "s#^ *DYNAMIC_VEG_OPTION *=.*# DYNAMIC_VEG_OPTION                = 8#" -e "s#^ *BTR_OPTION *=.*# BTR_OPTION                        = 4#" \
      -e "s#^ *SOIL_DATA_OPTION *=.*# SOIL_DATA_OPTION                  = 4#" \
      $N/src/hrldas/hrldas/run/examples/single_point/namelist.hrldas.single_point > namelist.hrldas
  sed -i "0,/^\/ *$/s##  SPATIAL_FILENAME = \"$R/forcing/soil_spatial.nc\"\n  PRESCRIBE_SOIL_MOISTURE = 1\n  LAIO_SSTAR = $LAIO_SSTAR\n  LAIO_SWILT = $LAIO_SWILT\n/#" namelist.hrldas
  for e in "$@"; do sed -i "$e" namelist.hrldas; done
  ./hrldas.exe > hrldas.log 2>&1; local rc=$?
  echo "== $C: exit $rc, DATE lines $(grep -c 'DATE=' hrldas.log), outputs $(ls out | grep -c LDASOUT)"
  [ $rc -ne 0 ] && { grep -v "k, \|Near) match\|^ *$" hrldas.log | tail -6; }
}
case_run full
case_run no_prescribe  "s#^ *PRESCRIBE_SOIL_MOISTURE *=.*#  PRESCRIBE_SOIL_MOISTURE = 0#"
case_run btr1          "s#^ *BTR_OPTION *=.*# BTR_OPTION = 1#"
case_run soil_table    "s#^ *SOIL_DATA_OPTION *=.*# SOIL_DATA_OPTION = 1#"
case_run dveg4         "s#^ *DYNAMIC_VEG_OPTION *=.*# DYNAMIC_VEG_OPTION = 4#"
case_run stock_opts    "s#^ *PRESCRIBE_SOIL_MOISTURE *=.*#  PRESCRIBE_SOIL_MOISTURE = 0#" "s#^ *BTR_OPTION *=.*# BTR_OPTION = 1#" \
                       "s#^ *SOIL_DATA_OPTION *=.*# SOIL_DATA_OPTION = 1#" "s#^ *DYNAMIC_VEG_OPTION *=.*# DYNAMIC_VEG_OPTION = 4#"
echo DEBUG_DONE
