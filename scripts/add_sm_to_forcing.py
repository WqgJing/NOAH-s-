"""TEST ONLY: add SMSURF / SMROOT to Bondville LDASIN forcing files, taken from a reference run's own output, so the
PRESCRIBE_SOIL_MOISTURE mechanics can be checked without L4 (1998 predates SMAP).
  SMSURF(t) = reference SOIL_M layer 1 (0-10 cm) at time t
  SMROOT(t) = reference SOIL_M thickness-weighted mean of layers 1-3 (0-100 cm, the L4 root-zone depth) at time t
Usage: python add_sm_to_forcing.py <forcing_dir> <reference_out_dir>
Writes the two variables into every <yyyymmddhhmm>.LDASIN_DOMAIN1 in place (dims Time, south_north, west_east)."""
import glob, os, sys
import numpy as np, pandas as pd, xarray as xr, netCDF4

fdir, rdir = sys.argv[1], sys.argv[2]
parts = []
for f in sorted(glob.glob(f"{rdir}/*.LDASOUT_DOMAIN1")):
    with xr.open_dataset(f, decode_times=False) as ds:
        t = pd.to_datetime([s.decode() if isinstance(s, bytes) else str(s) for s in ds["Times"].values], format="%Y-%m-%d_%H:%M:%S")
        sm = ds["SOIL_M"].values[:, 0, :, 0]                      # (Time, soil layers)
        parts.append(pd.DataFrame({"surf": sm[:, 0], "root": (0.1 * sm[:, 0] + 0.3 * sm[:, 1] + 0.6 * sm[:, 2]) / 1.0}, index=t))
ref = pd.concat(parts).sort_index(); ref = ref[~ref.index.duplicated()]
files = sorted(glob.glob(f"{fdir}/*.LDASIN_DOMAIN1"))
nmiss = 0
for f in files:
    t = pd.to_datetime(os.path.basename(f)[:12], format="%Y%m%d%H%M")
    if t in ref.index:
        s, r = float(ref.at[t, "surf"]), float(ref.at[t, "root"])
    else:                                                          # first forcing time (no output yet): nearest later value
        j = ref.index.searchsorted(t); s, r = float(ref["surf"].iloc[j]), float(ref["root"].iloc[j]); nmiss += 1
    with netCDF4.Dataset(f, "a") as nc:
        for name, val in (("SMSURF", s), ("SMROOT", r)):
            v = nc.variables[name] if name in nc.variables else nc.createVariable(name, "f4", ("Time", "south_north", "west_east"))
            v[0, 0, 0] = val; v.units = "m3/m3"
print(f"added SMSURF/SMROOT to {len(files)} files ({nmiss} filled from the next output time); "
      f"SMSURF {ref.surf.min():.3f}-{ref.surf.max():.3f}, SMROOT {ref.root.min():.3f}-{ref.root.max():.3f}")
