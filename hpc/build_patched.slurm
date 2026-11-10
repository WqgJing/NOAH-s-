#!/bin/bash
#SBATCH --job-name=hrldas_patched
#SBATCH --output=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah/logs/build_patched_%j.out
#SBATCH --account=gts-rbras6-coda20 --partition=cpu-small -q inferno
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:45:00
# Build the PATCHED HRLDAS/Noah-MP in build/hrldas (src/ + patches/noahmp_sstar.patch). Same compiler, flags and
# netCDF paths as the stock build (hpc/build_hrldas.slurm), so differences in results come only from the patch.
set -e
N=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah
SRC=$N/build/hrldas/hrldas
source /etc/profile.d/z00_lmod_pace.sh   # lmod is not initialised in non-login batch shells
source ~/.bashrc
module load gcc/12.3.0 mvapich2/2.3.7-1 netcdf-c/4.9.2-mva2-hdf5-1.14 netcdf-fortran/4.6.1-mva2-hdf5-1.14 || true
gfortran --version | head -1
NFROOT=$(nf-config --prefix); NCLIBS=$(nc-config --libs)
cd $SRC
echo 6 | ./configure > /dev/null                                  # 6 = Linux gfortran/gcc serial (as stock)
sed -i "s#^ *NETCDFMOD *=.*# NETCDFMOD      =      -I$NFROOT/include#" user_build_options
sed -i "s#^ *NETCDFLIB *=.*# NETCDFLIB      =      -L$NFROOT/lib -lnetcdff $NCLIBS#" user_build_options
diff <(grep -v "^#" $N/src/hrldas/hrldas/user_build_options) <(grep -v "^#" user_build_options) && echo "build options identical to stock"
make > $N/logs/build_patched_make_${SLURM_JOB_ID}.log 2>&1 || echo "WARNING: top-level make returned non-zero (see make log)"
grep -i -E "error" $N/logs/build_patched_make_${SLURM_JOB_ID}.log | grep -v -i "errorhandle\|error_handle\|ErrorHandleMod" | head -20 || true
ls -la run/hrldas.exe
echo BUILD_OK
