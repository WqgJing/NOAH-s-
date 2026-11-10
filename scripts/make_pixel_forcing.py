"""HRLDAS inputs for one Noah-MP test pixel (noah/data/l4_points.npz: BV, VR, KG) from our own data.
Usage: python make_pixel_forcing.py <NAME> <YEAR0> <YEAR1> <OUTDIR>
Writes into OUTDIR: hourly <yyyymmddhh>.LDASIN_DOMAIN1 from YEAR0-01-01 00 to (YEAR1+1)-01-01 00 UTC, hrldas_setup.nc,
soil_spatial.nc, pixel_info.txt (values the run script puts in the namelist: LAIO_SSTAR, LAIO_SWILT, KDAY).

Forcing (ERA5 hourly, clima/forcing_global/era5_hourly/<year>/<var>.npy, column = the pixel's index in pixels.npz):
  T2D = t2m [K]; PSFC = sp [Pa]; Q2D from d2m and sp [kg/kg]; U2D, V2D = u10, v10 [m/s] (ZLVL = 10 m, as the stock example)
  SWDOWN, LWDOWN [W/m2]: ERA5 mean fluxes are means over the hour ENDING at the stamp; the value at t is centred on t:
                          0.5 * (m(t) + m(t+1h)). HRLDAS interpolates SW/LW linearly between forcing times.
  RAINRATE [mm/s] at t = mtpr(t+1h): HRLDAS holds precipitation of the file at t over (t, t+1h].
  LAI  = Yuan v6.1 8-day (pixels/lai_yuan_v061.npz), composite centred at doy+3.5, linear in time   (DYNAMIC_VEG_OPTION = 8)
  SMROOT = SMAP L4 root zone from the LOCAL daily files (sbi_train/runs/l4/l4_rootzone_<year>.npy; the map's product),
           daily values at 12:00 UTC, linear in time                                      (PRESCRIBE_SOIL_MOISTURE = 1)
  SMSURF = SMAP L4 (SPL4SMGP/008) sm_surface, 3-hourly from GEE (no local L4 surface), linear in time
Soil (SOIL_DATA_OPTION = 4, soil_spatial.nc), all 4 layers the same:
  smcmax = L4 porosity, bexp = L4 b, psisat = |L4 psi_s| (Campbell, as L4 / the bank),
  smcwlt = smcdry = theta(-153 m) and smcref = theta(-3.36 m) on that Campbell curve,
  dksat, quartz and the runoff / VIC / irrigation / tile placeholders from the Noah-MP table for ISLTYP (loam, 6; theta is prescribed,
  so drainage parameters do not act; quartz sets soil thermal conductivity), dwsat = b * dksat * psisat / smcmax.
Laio parameters: s* = bank theta_star_ref, s_w = bank theta_w (the L4 reference values used elsewhere in the project)."""
import os, re, sys
import numpy as np, pandas as pd, netCDF4

NAME, Y0, Y1, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
D = "/storage/project/r-rbras6-0/jjason3/proposal_drought"
ERA = f"{D}/clima/forcing_global/era5_hourly"
TBL = f"{D}/noah/src/hrldas/hrldas/run/NoahmpTable.TBL"
ISLTYP = 6                                                       # loam (table placeholders only, see docstring)
DZS = np.array([0.10, 0.30, 0.60, 1.00]); ZS = np.array([0.05, 0.25, 0.70, 1.50])
os.makedirs(OUT, exist_ok=True)

L = np.load(f"{D}/noah/data/l4_points.npz", allow_pickle=True)
j = list(L["name"]).index(NAME)
pix, bank, lat, lon = int(L["pix"][j]), int(L["bank"][j]), float(L["lat"][j]), float(L["lon"][j])
igbp, poros, bee, psis = int(L["igbp"][j]), float(L["poros"][j]), float(L["bee"][j]), abs(float(L["psis"][j]))
B = np.load(f"{D}/sbi_train/runs/vegsbi/bank.npz")
sstar, swilt_bank = float(B["theta_star_ref"][bank]), float(B["theta_w"][bank])
camp = lambda psi_m: min(poros, poros * (psi_m / psis) ** (-1.0 / bee))   # Campbell theta at suction psi_m [m]
smcwlt, smcref = camp(153.0), camp(3.36)
print(f"{NAME}: pixel {pix} bank {bank} ({lat:.3f}, {lon:.3f}) IGBP {igbp}; poros {poros:.3f} b {bee:.2f} psis {psis:.3f} m; "
      f"theta(-153 m) {smcwlt:.3f} (bank theta_w {swilt_bank:.3f}); theta(-3.36 m) {smcref:.3f}; s* {sstar:.3f}", flush=True)

# ---- hourly time axis (UTC) and ERA5
t_out = pd.date_range(f"{Y0}-01-01 00:00", f"{Y1 + 1}-01-01 00:00", freq="h")
t_era = pd.date_range(f"{Y0}-01-01 00:00", f"{Y1 + 1}-12-31 23:00", freq="h")
V = {}
for v in ("t2m", "d2m", "sp", "u10", "v10", "msdwswrf", "msdwlwrf", "mtpr"):
    parts = [np.load(f"{ERA}/{y}/{v}.npy", mmap_mode="r")[:, pix].astype(np.float64) for y in range(Y0, Y1 + 2)]
    V[v] = pd.Series(np.concatenate(parts), index=t_era).reindex(t_out.union([t_out[-1] + pd.Timedelta("1h")]))
for v, s in V.items():
    if s.isna().any(): raise RuntimeError(f"ERA5 {v}: {s.isna().sum()} missing hours")
td, ps = V["d2m"].to_numpy() - 273.15, V["sp"].to_numpy()
e = 611.2 * np.exp(17.67 * td / (td + 243.5))                    # vapour pressure from dewpoint [Pa] (Bolton 1980)
q = 0.622 * e / (ps - 0.378 * e)
centre = lambda m: 0.5 * (m[:-1] + m[1:])                         # value at t = mean over (t-30 min, t+30 min]
F = {"T2D": V["t2m"].to_numpy()[:-1], "Q2D": q[:-1], "PSFC": ps[:-1], "U2D": V["u10"].to_numpy()[:-1], "V2D": V["v10"].to_numpy()[:-1],
     "SWDOWN": np.maximum(0.0, centre(V["msdwswrf"].to_numpy())), "LWDOWN": centre(V["msdwlwrf"].to_numpy()),
     "RAINRATE": np.maximum(0.0, V["mtpr"].to_numpy()[1:])}

# ---- LAI (8-day, centred) and L4 (3-hourly) -> hourly
Z = np.load(f"{D}/clima/forcing_global/pixels/lai_yuan_v061.npz"); doy = Z["doy"]
lt, lv = [], []
for y in range(Y0 - 1, Y1 + 2):
    if f"leaf_area_index_{y}" not in Z.files: continue
    lt += list(pd.Timestamp(f"{y}-01-01") + pd.to_timedelta(doy - 1 + 3.5, "D")); lv += list(Z[f"leaf_area_index_{y}"][pix])
lai = pd.Series(lv, index=pd.DatetimeIndex(lt)).sort_index().dropna()
tsec = lambda ti: ti.values.astype("datetime64[s]").astype(np.int64)
F["LAI"] = np.maximum(0.0, np.interp(tsec(t_out), tsec(lai.index), lai.to_numpy()))
tl = pd.DatetimeIndex(L["time"]); ok = np.isfinite(L["sm_surface"][:, j])
F["SMSURF"] = np.interp(tsec(t_out), tsec(tl[ok]), L["sm_surface"][ok, j])
print(f"L4 surface (GEE 008, 3-hourly): first {tl[0]}, hours of day {sorted(set(tl[:16].hour))}; missing 3-h steps {(~ok).sum()}", flush=True)
# Root zone from the LOCAL daily L4 files, i.e. the product the Laio map uses (identical to GEE 008 daily means 2018-2022,
# within 3e-4 in 2023; 2024 local = version 007). Daily means of UTC days, placed at 12:00 UTC, linear in time, gaps bridged.
grid = pd.read_csv("/storage/project/r-rbras6-0/jjason3/forward_sensitivity_lib/data/global/global_grid_025deg.csv").sort_values("grid_id")
grow = int(np.searchsorted(grid.grid_id.to_numpy(), int(B["grid_id"][bank])))
assert int(grid.grid_id.iloc[grow]) == int(B["grid_id"][bank])
rt, rv = [], []
for y in range(Y0 - 1, Y1 + 2):
    fn = f"{D}/sbi_train/runs/l4/l4_rootzone_{y}.npy"
    if not os.path.exists(fn): continue
    a = np.asarray(np.load(fn, mmap_mode="r")[grow], float)
    rt += list(pd.date_range(f"{y}-01-01 12:00", periods=len(a), freq="D")); rv += list(a)
root = pd.Series(rv, index=pd.DatetimeIndex(rt)).dropna()
F["SMROOT"] = np.interp(tsec(t_out), tsec(root.index), root.to_numpy())
nmiss = int(((root.index >= t_out[0]) & (root.index <= t_out[-1])).sum())
print(f"L4 root zone (local daily, grid row {grow}): {nmiss} valid days in the run period of {(t_out[-1] - t_out[0]).days}", flush=True)
for k, a in F.items():
    print(f"  {k:9s} min {np.min(a):10.4g} mean {np.mean(a):10.4g} max {np.max(a):10.4g}", flush=True)

# ---- Noah-MP table values for the placeholder soil class
txt = open(TBL).read()
def block(header):                                               # namelist block ends at a line that is only "/" (comments contain "/")
    b = txt[txt.index(header):]; return b[:re.search(r"^\s*/\s*$", b, re.M).start()]
blk = block("&noahmp_soil_stas_parameters")
def tab(name):
    m = re.search(rf"^\s*{name}\s*=\s*([^\n!]+)", blk, re.M); return float(m.group(1).split(",")[ISLTYP - 1])
dksat, quartz = tab("SATDK"), tab("QTZ"); dwsat = bee * dksat * psis / poros
gen = block("&noahmp_general_parameters")
gval = lambda name, dflt: float(re.search(rf"^\s*{name}\s*=\s*([^\n!,]+)", gen, re.M).group(1)) if re.search(rf"^\s*{name}\s*=", gen, re.M) else dflt

def dims2(nc):
    nc.createDimension("Time", None); nc.createDimension("south_north", 1); nc.createDimension("west_east", 1)

# ---- spatial soil file
with netCDF4.Dataset(f"{OUT}/soil_spatial.nc", "w") as nc:
    nc.createDimension("soil_layers_stag", 4); nc.createDimension("south_north", 1); nc.createDimension("west_east", 1)
    for name, val in (("bexp", bee), ("smcdry", smcwlt), ("smcwlt", smcwlt), ("smcref", smcref), ("smcmax", poros),
                      ("dksat", dksat), ("dwsat", dwsat), ("psisat", psis), ("quartz", quartz)):
        nc.createVariable(name, "f4", ("soil_layers_stag", "south_north", "west_east"))[:] = val
    two = {"refdk": gval("REFDK_DATA", 2.0e-6), "refkdt": gval("REFKDT_DATA", 3.0)}
    for name in ("BVIC", "AXAJ", "BXAJ", "XXAJ", "BDVIC", "GDVIC", "BBVIC"):
        two[name] = tab(name)
    for name in ("IRR_FRAC", "IRR_HAR", "IRR_LAI", "IRR_MAD", "FILOSS", "SPRIR_RATE", "MICIR_RATE", "FIRTFAC", "IR_RAIN",
                 "KLAT_FAC", "TDSMC_FAC", "TD_DC", "TD_DCOEF", "TD_DDRAIN", "TD_RADI", "TD_SPAC"):
        two[name] = 0.0                                          # irrigation / tile drainage are off; read but not used
    for name, val in two.items():
        nc.createVariable(name, "f4", ("south_north", "west_east"))[:] = val

# ---- setup file (as create_point_data.f90)
fveg = lambda l: 100.0 * (1.0 - np.exp(-0.52 * (l + max(0.05, 0.1 * l))))
with netCDF4.Dataset(f"{OUT}/hrldas_setup.nc", "w") as nc:
    dims2(nc); nc.createDimension("soil_layers_stag", 4)
    def put(n, val, dt="f4"):                                    # write record 0 explicitly: Time is unlimited (length 0),
        nc.createVariable(n, dt, ("Time", "south_north", "west_east"))[0, 0, 0] = val   # so var[:] = x would write nothing
    put("XLAT", lat); put("XLONG", lon); put("TMN", float(np.mean(F["T2D"]))); put("HGT", 0.0); put("SEAICE", 0.0)
    put("SHDMAX", fveg(F["LAI"].max())); put("SHDMIN", fveg(F["LAI"].min())); put("LAI", F["LAI"][0])
    put("XLAND", 1, "i4"); put("IVGTYP", igbp, "i4"); put("ISLTYP", ISLTYP, "i4")
    put("SNOW", 0.0); put("SNODEP", 0.0); put("CANWAT", 0.0); put("TSK", F["T2D"][0])
    nc.createVariable("DZS", "f4", ("Time", "soil_layers_stag"))[0, :] = DZS
    nc.createVariable("ZS", "f4", ("Time", "soil_layers_stag"))[0, :] = ZS
    nc.createVariable("TSLB", "f4", ("Time", "soil_layers_stag", "south_north", "west_east"))[0, :, 0, 0] = \
        np.linspace(F["T2D"][0], float(np.mean(F["T2D"])), 4)
    nc.createVariable("SMOIS", "f4", ("Time", "soil_layers_stag", "south_north", "west_east"))[0, :, 0, 0] = \
        [F["SMSURF"][0], F["SMROOT"][0], F["SMROOT"][0], F["SMROOT"][0]]
    for k, v in (("GRID_ID", 1), ("ISWATER", 17), ("ISURBAN", 13), ("ISICE", 15), ("MMINLU", "MODIFIED_IGBP_MODIS_NOAH"),
                 ("DX", 0.0), ("DY", 0.0), ("TRUELAT1", 0.0), ("TRUELAT2", 0.0), ("STAND_LON", 0.0), ("MAP_PROJ", 0),
                 ("TITLE", "noah/scripts/make_pixel_forcing.py")):
        nc.setncattr(k, v)

# ---- the same inputs as arrays (kept after the run; the Python port reads these instead of the LDASIN files)
np.savez(f"{OUT}/forcing_arrays.npz", time=t_out.values.astype("datetime64[s]"), **{k: np.asarray(a, np.float64) for k, a in F.items()},
         lat=lat, lon=lon, igbp=igbp, isltyp=ISLTYP, tmn=float(np.mean(F["T2D"])), dzs=DZS, zs=ZS,
         tslb0=np.linspace(F["T2D"][0], float(np.mean(F["T2D"])), 4), smois0=[F["SMSURF"][0]] + [F["SMROOT"][0]] * 3,
         bexp=bee, smcmax=poros, psisat=psis, smcwlt=smcwlt, smcref=smcref, dksat=dksat, dwsat=dwsat, quartz=quartz,
         laio_sstar=sstar, laio_swilt=swilt_bank)

# ---- hourly forcing files
for i, t in enumerate(t_out):
    with netCDF4.Dataset(f"{OUT}/{t:%Y%m%d%H}.LDASIN_DOMAIN1", "w") as nc:
        dims2(nc)
        for k, a in F.items():
            nc.createVariable(k, "f4", ("Time", "south_north", "west_east"))[0, 0, 0] = a[i]
        nc.setncattr("TITLE", "noah/scripts/make_pixel_forcing.py")
with netCDF4.Dataset(f"{OUT}/hrldas_setup.nc") as nc:            # read back: every setup variable must hold a real value
    bad = [k for k, v in nc.variables.items() if np.ma.is_masked(v[:]) or v[:].size == 0]
    print(f"setup check: XLAT {float(nc['XLAT'][0, 0, 0]):.3f}, IVGTYP {int(nc['IVGTYP'][0, 0, 0])}, ISLTYP {int(nc['ISLTYP'][0, 0, 0])}, "
          f"TMN {float(nc['TMN'][0, 0, 0]):.1f}; empty variables: {bad}", flush=True)
    if bad: raise RuntimeError(f"setup file has empty variables: {bad}")
kday = (t_out[-1] - t_out[0]).days
with open(f"{OUT}/pixel_info.txt", "w") as f:
    f.write(f"NAME={NAME}\nLAT={lat}\nLON={lon}\nIGBP={igbp}\nLAIO_SSTAR={sstar:.4f}\nLAIO_SWILT={swilt_bank:.4f}\nKDAY={kday}\n"
            f"START={t_out[0]:%Y-%m-%d}\n")
print(f"wrote {len(t_out)} forcing files ({t_out[0]} .. {t_out[-1]}), setup, spatial soil; KDAY {kday}", flush=True)
