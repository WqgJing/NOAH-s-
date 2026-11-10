"""Set up a variant of a Fortran reference run (runs/pyref_<NAME>) with other Laio / nuisance parameters, to test
code paths the stock reference does not reach (e.g. BTRAN < 1). Used by hpc/pyref_variant.slurm.

Usage: python make_pyref_variant.py <SRC run dir> <DST run dir> <LAIO_SSTAR> <LAIO_SWILT> [GS_SCALE RSURF_SCALE Z0_SCALE]
Copies forcing/{hrldas_setup.nc, soil_spatial.nc, forcing_arrays.npz} and *.TBL from SRC, rewrites the hourly
LDASIN files of the run window (start .. start + KDAY + 1 h) from forcing_arrays.npz exactly as
make_pixel_forcing.py writes them (float32), and writes DST/namelist.hrldas = SRC namelist with the new paths and
parameters. Inputs are otherwise identical to SRC.
"""
import glob
import os
import re
import shutil
import sys

import netCDF4
import numpy as np
import pandas as pd

src, dst, sstar, swilt = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
scales = [float(x) for x in sys.argv[5:8]] if len(sys.argv) > 5 else None
os.makedirs(f"{dst}/forcing", exist_ok=True)
os.makedirs(f"{dst}/out", exist_ok=True)
for f in ("hrldas_setup.nc", "soil_spatial.nc", "forcing_arrays.npz"):
    shutil.copy(f"{src}/forcing/{f}", f"{dst}/forcing/{f}")
for f in glob.glob(f"{src}/*.TBL"):
    shutil.copy(f, dst)

nl = open(f"{src}/namelist.hrldas").read()
get = lambda k: int(re.search(rf"^\s*{k}\s*=\s*(\d+)", nl, re.M | re.I).group(1))
start = pd.Timestamp(get("START_YEAR"), get("START_MONTH"), get("START_DAY"), get("START_HOUR"))
end = start + pd.Timedelta(days=get("KDAY")) + pd.Timedelta(hours=1)
for key, val in (("HRLDAS_SETUP_FILE", f'"{dst}/forcing/hrldas_setup.nc"'), ("INDIR", f'"{dst}/forcing/"'),
                 ("OUTDIR", f'"{dst}/out/"'), ("SPATIAL_FILENAME", f'"{dst}/forcing/soil_spatial.nc"'),
                 ("LAIO_SSTAR", f"{sstar:.4f}"), ("LAIO_SWILT", f"{swilt:.4f}")):
    nl, n = re.subn(rf"^(\s*{key}\s*=).*$", rf"\g<1> {val}", nl, flags=re.M | re.I)
    if n != 1:
        raise RuntimeError(f"{key}: {n} matches in the source namelist")
if scales:
    add = "".join(f"  {k} = {s}\n" for k, s in zip(("NUIS_GS_SCALE", "NUIS_RSURF_SCALE", "NUIS_Z0_SCALE"), scales))
    nl = re.sub(r"^[ \t]*/[ \t]*$", add + "/", nl, count=1, flags=re.M)   # keep the newline after "/"
open(f"{dst}/namelist.hrldas", "w").write(nl)

Z = np.load(f"{dst}/forcing/forcing_arrays.npz")
t = pd.DatetimeIndex(Z["time"])
keys = ("T2D", "Q2D", "PSFC", "U2D", "V2D", "SWDOWN", "LWDOWN", "RAINRATE", "LAI", "SMSURF", "SMROOT")
sel = np.flatnonzero((t >= start) & (t <= end))
for i in sel:
    with netCDF4.Dataset(f"{dst}/forcing/{t[i]:%Y%m%d%H}.LDASIN_DOMAIN1", "w") as nc:
        nc.createDimension("Time", None); nc.createDimension("south_north", 1); nc.createDimension("west_east", 1)
        for k in keys:
            nc.createVariable(k, "f4", ("Time", "south_north", "west_east"))[0, 0, 0] = Z[k][i]
        nc.setncattr("TITLE", "noah/scripts/make_pixel_forcing.py")
print(f"variant {dst}: s* {sstar} s_w {swilt} scales {scales}; {len(sel)} LDASIN files {t[sel[0]]} .. {t[sel[-1]]}")
