#!/bin/bash
# Recreate the patched build copy from pristine NCAR code: build/hrldas = src/hrldas + patches/noahmp_sstar.patch
# (src/ is never edited; run after a fresh clone with submodules, then sbatch hpc/build_patched.slurm)
set -e
N=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah
cd $N
rm -rf build/hrldas; mkdir -p build
rsync -a --exclude '.git' --exclude '*.o' --exclude '*.mod' --exclude '*.exe' --exclude '*.a' --exclude 'user_build_options' src/hrldas/ build/hrldas/
# drop preprocessed copies (.f90 with a .F90 twin); the Makefiles regenerate them
for f in $(find build/hrldas -name '*.f90'); do [ -f "${f%.f90}.F90" ] && rm -f "$f"; done
patch -p2 -d build/hrldas < patches/noahmp_sstar.patch
echo "build/hrldas ready (patched)"
