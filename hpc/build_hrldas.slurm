#!/bin/bash
#SBATCH --job-name=hrldas_build
#SBATCH --output=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah/logs/build_%j.out
#SBATCH --account=gts-rbras6-coda20 --partition=cpu-small -q inferno
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:45:00
# Build stock HRLDAS v5.2.1 + Noah-MP (serial gfortran, configure option 6) and the single-point forcing tool.
# Only change from upstream: netCDF include/lib paths in user_build_options (PACE modules instead of NCAR /glade paths).
set -e
N=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah
SRC=$N/src/hrldas/hrldas
source /etc/profile.d/z00_lmod_pace.sh   # lmod is not initialised in non-login batch shells
source ~/.bashrc
module load gcc/12.3.0 mvapich2/2.3.7-1 netcdf-c/4.9.2-mva2-hdf5-1.14 netcdf-fortran/4.6.1-mva2-hdf5-1.14 || true
module list 2>&1 | tail -8
gfortran --version | head -1
NFROOT=$(nf-config --prefix); NCLIBS=$(nc-config --libs)
echo "netcdf-fortran: $NFROOT"; echo "nc-config --libs: $NCLIBS"

cd $SRC
echo 6 | ./configure > /dev/null                                  # 6 = Linux gfortran/gcc serial
sed -i "s#^ *NETCDFMOD *=.*# NETCDFMOD      =      -I$NFROOT/include#" user_build_options
sed -i "s#^ *NETCDFLIB *=.*# NETCDFLIB      =      -L$NFROOT/lib -lnetcdff $NCLIBS#" user_build_options
grep -E "COMPILERF90|F90FLAGS|NETCDFMOD|NETCDFLIB" user_build_options

# Top-level make builds hrldas.exe before HRLDAS_forcing (GRIB tool may need jasper; we do not use it).
make > $N/logs/build_make_${SLURM_JOB_ID}.log 2>&1 || echo "WARNING: top-level make returned non-zero (see build_make log)"
ls -la run/hrldas.exe

cd $SRC/HRLDAS_forcing/run/examples/single_point
make > $N/logs/build_point_${SLURM_JOB_ID}.log 2>&1
ls -la create_point_data.exe
echo BUILD_OK
