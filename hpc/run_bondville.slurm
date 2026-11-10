#!/bin/bash
#SBATCH --job-name=bondville
#SBATCH --output=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah/logs/bondville_%x_%j.out
#SBATCH --account=gts-rbras6-coda20 --partition=cpu-small -q inferno
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=01:00:00
# Reproduce the stock HRLDAS single-point example (Bondville, IL, 1998; 30-min forcing).
#   CONFIG=shipped  : run/examples/single_point/namelist.hrldas.single_point as distributed
#   CONFIG=tutorial : the options shown in tutorial/Note1_Single_Point.ipynb screenshots
#                     (SNOW_ALBEDO 2, TEMP_TIME_SCHEME 1, SURFACE_RESISTANCE 1, CROP 1)
# Both: start 1998-01-02 00:00 (the forcing starts then; the tutorial also used this), KDAY 364 (forcing ends 1999-01-01 00:00),
# paths moved to noah/runs/, SPLIT_OUTPUT_COUNT 48 (daily files; does not change results).
set -e
CONFIG=${CONFIG:-shipped}
N=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah
SRC=$N/src/hrldas/hrldas
EX=$SRC/HRLDAS_forcing/run/examples/single_point
RUN=$N/runs/bondville_$CONFIG
source /etc/profile.d/z00_lmod_pace.sh   # lmod is not initialised in non-login batch shells
source ~/.bashrc
module load gcc/12.3.0 mvapich2/2.3.7-1 netcdf-c/4.9.2-mva2-hdf5-1.14 netcdf-fortran/4.6.1-mva2-hdf5-1.14 || true

rm -rf $RUN; mkdir -p $RUN/forcing $RUN/out
cd $RUN/forcing
cp $EX/bondville.dat .
SECONDS=0; $EX/create_point_data.exe > $RUN/create_point_data.log 2>&1
echo "create_point_data: ${SECONDS}s, LDASIN files: $(ls | grep -c LDASIN_DOMAIN1)"

cd $RUN
cp $SRC/run/*.TBL .
SED=(-e "s#^ *HRLDAS_SETUP_FILE *=.*# HRLDAS_SETUP_FILE = \"$RUN/forcing/hrldas_setup_single_point.nc\"#"
     -e "s#^ *INDIR *=.*# INDIR = \"$RUN/forcing/\"#"
     -e "s#^ *OUTDIR *=.*# OUTDIR = \"$RUN/out/\"#"
     -e "s#^ *START_DAY *=.*# START_DAY   = 02#"
     -e "s#^ *KDAY *=.*# KDAY = 364#"
     -e "s#^ *SPLIT_OUTPUT_COUNT *=.*# SPLIT_OUTPUT_COUNT = 48#")
if [ "$CONFIG" = "tutorial" ]; then
  SED+=(-e "s#^ *SNOW_ALBEDO_OPTION *=.*# SNOW_ALBEDO_OPTION                = 2#"
        -e "s#^ *TEMP_TIME_SCHEME_OPTION *=.*# TEMP_TIME_SCHEME_OPTION           = 1#"
        -e "s#^ *SURFACE_RESISTANCE_OPTION *=.*# SURFACE_RESISTANCE_OPTION         = 1#"
        -e "s#^ *CROP_OPTION *=.*# CROP_OPTION                       = 1#")
fi
sed "${SED[@]}" $SRC/run/examples/single_point/namelist.hrldas.single_point > namelist.hrldas
echo "== namelist changes vs shipped ($CONFIG)"; diff $SRC/run/examples/single_point/namelist.hrldas.single_point namelist.hrldas || true
ln -sf $SRC/run/hrldas.exe .
SECONDS=0; ./hrldas.exe > hrldas.log 2>&1 || { echo "hrldas.exe FAILED"; tail -30 hrldas.log; exit 1; }
echo "hrldas.exe: ${SECONDS}s for 1 point-year"
grep -m 16 '\*\*\*DATE=' hrldas.log | cut -c1-60
echo "output files: $(ls out | wc -l)"
# keep the run directory small: drop the 17,520 forcing files (re-creatable from bondville.dat in seconds)
rm -f $RUN/forcing/*.LDASIN_DOMAIN1
echo RUN_OK
