#!/bin/bash
# Export our source changes (build/hrldas vs pristine src/hrldas) as one patch: patches/noahmp_sstar.patch
# Only source files are compared (build products, generated .f90 and local build options are excluded).
set -e
N=/storage/project/r-rbras6-0/jjason3/proposal_drought/noah
cd $N
diff -ruN -x '.git' -x '*.o' -x '*.mod' -x '*.f90' -x '*.exe' -x '*.a' -x 'user_build_options' \
     -x '*.log' -x '*.nc' src/hrldas build/hrldas > patches/noahmp_sstar.patch || true   # diff exits 1 when files differ
echo "files changed:"; grep '^+++ ' patches/noahmp_sstar.patch | awk '{print "  " $2}'
echo "lines: +$(grep -c '^+[^+]' patches/noahmp_sstar.patch) -$(grep -c '^-[^-]' patches/noahmp_sstar.patch)"
