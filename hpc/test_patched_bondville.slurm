#!/bin/bash
#SBATCH --job-name=test_patched
#SBATCH --output=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah/logs/test_patched_%j.out
#SBATCH --account=gts-rbras6-coda20 --partition=cpu-small -q inferno
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=01:30:00
# Tests of the patched Noah-MP (build/hrldas) on NCAR's Bondville 1998 example, shipped options, same run settings as
# hpc/run_bondville.slurm (start 1998-01-02, KDAY 364, daily output). Cases:
#   R  : all switches off                         -> must equal the stock run runs/bondville_shipped bit for bit
#   B4 : BTR_OPTION 4 (Laio), s* = REFSMC 0.387, s_w = WLTSMC 0.120 (soil class 8 table values)
#   P  : B4 + PRESCRIBE_SOIL_MOISTURE 1, SMSURF/SMROOT taken from the stock run's own SOIL_M (test only)
#   S1/S2/S3 : NUIS_GS_SCALE / NUIS_RSURF_SCALE / NUIS_Z0_SCALE = 2.0 (direction checks)
set -e
N=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah
STOCK=$N/src/hrldas/hrldas
PATCHED=$N/build/hrldas/hrldas
EX=$STOCK/HRLDAS_forcing/run/examples/single_point
T=$N/runs/test_patched
source /etc/profile.d/z00_lmod_pace.sh   # lmod is not initialised in non-login batch shells
source ~/.bashrc
module load gcc/12.3.0 mvapich2/2.3.7-1 netcdf-c/4.9.2-mva2-hdf5-1.14 netcdf-fortran/4.6.1-mva2-hdf5-1.14 || true
conda activate gw_viz

rm -rf $T; mkdir -p $T/forcing
cd $T/forcing; cp $EX/bondville.dat .; $EX/create_point_data.exe > /dev/null 2>&1
echo "forcing files: $(ls | grep -c LDASIN_DOMAIN1)"

base_namelist () {   # $1 = output dir; stock example namelist with the same run settings as the stock Bondville run
  sed -e "s#^ *HRLDAS_SETUP_FILE *=.*# HRLDAS_SETUP_FILE = \"$T/forcing/hrldas_setup_single_point.nc\"#" \
      -e "s#^ *INDIR *=.*# INDIR = \"$T/forcing/\"#" -e "s#^ *OUTDIR *=.*# OUTDIR = \"$1/\"#" \
      -e "s#^ *START_DAY *=.*# START_DAY   = 02#" -e "s#^ *KDAY *=.*# KDAY = 364#" \
      -e "s#^ *SPLIT_OUTPUT_COUNT *=.*# SPLIT_OUTPUT_COUNT = 48#" \
      $STOCK/run/examples/single_point/namelist.hrldas.single_point
}
run_case () {   # $1 = case name, $2.. = extra namelist lines (inserted before the closing "/")
  local C=$1; shift; local D=$T/$C; mkdir -p $D/out; cd $D
  cp $PATCHED/run/*.TBL .; ln -sf $PATCHED/run/hrldas.exe .
  base_namelist $D/out > namelist.tmp
  if [ $# -gt 0 ]; then
    for L in "$@"; do
      key=$(echo "$L" | cut -d= -f1 | tr -d ' ')
      if grep -q -i "^ *$key *=" namelist.tmp; then sed -i "s#^ *$key *=.*# $L#I" namelist.tmp; else sed -i "0,/^\/ *$/s##  $L\n/#" namelist.tmp; fi
    done
  fi
  mv namelist.tmp namelist.hrldas
  SECONDS=0; ./hrldas.exe > hrldas.log 2>&1 || { echo "$C: hrldas.exe FAILED"; tail -20 hrldas.log; return 1; }
  echo "$C: ${SECONDS}s; $(grep 'SFLXtot' hrldas.log | tail -1 | sed 's/  */ /g')"
}

run_case R
run_case B4 "BTR_OPTION = 4" "LAIO_SSTAR = 0.387" "LAIO_SWILT = 0.120"
python -u $N/scripts/add_sm_to_forcing.py $T/forcing $N/runs/bondville_shipped/out
run_case P  "BTR_OPTION = 4" "LAIO_SSTAR = 0.387" "LAIO_SWILT = 0.120" "PRESCRIBE_SOIL_MOISTURE = 1"
run_case S1 "NUIS_GS_SCALE = 2.0"
run_case S2 "NUIS_RSURF_SCALE = 2.0"
run_case S3 "NUIS_Z0_SCALE = 2.0"
echo "== namelist differences vs R:"; for C in B4 P S1 S2 S3; do echo "-- $C"; diff $T/R/namelist.hrldas $T/$C/namelist.hrldas | grep "^>" || true; done
rm -f $T/forcing/*.LDASIN_DOMAIN1
cd $N; python -u scripts/test_patched_compare.py
echo TEST_DONE
