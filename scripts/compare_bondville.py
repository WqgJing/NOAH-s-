"""Did we reproduce NCAR's single-point example (Bondville, IL, 1998)?

References (from src/hrldas/tutorial/Note1_Single_Point.ipynb, extracted to runs/reference/):
  1. run log screenshot: per-step printout "***DATE=<date> TSLB(layer 1) LAI" for 1998-01-02 00:00-07:30 (exact numbers).
  2. ncview screenshot of LH in file 199806200030.LDASOUT_DOMAIN1: 48 half-hourly values from 1998-06-20 00:30 UTC,
     digitised by eye from the PNG (about +-10 W m-2).
The tutorial was made with an older Noah-MP release than ours (v5.2.1; its namelist lacks several v5.2 options), so small
differences are expected. Two runs: CONFIG=shipped (namelist as distributed) and CONFIG=tutorial (options in the screenshots).
Out: printed comparison, runs/figs/bondville_reproduce.png (line plots only)."""
import glob, re
import numpy as np, pandas as pd, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

N = "/storage/project/r-rbras6-0/jjason3/proposal_drought/noah"
CONFIGS = ["shipped", "tutorial"]

# reference 1: tutorial run log, 1998-01-02 00:00 to 07:30 (TSLB layer 1 [K])
REF_T = pd.date_range("1998-01-02 00:00", periods=16, freq="30min")
REF_TSLB = np.array([266.10001, 266.80774, 267.76233, 268.69067, 269.49200, 270.14383, 270.66016, 271.02570,
                     271.28119, 271.46808, 271.57654, 271.71628, 271.87170, 272.03931, 272.19501, 272.33582])
# reference 2: tutorial ncview LH [W m-2], 1998-06-20 00:30 UTC + 30 min x i (digitised, +-10)
REF_LH = np.array([57, 8, 5, 5, 5, 5, 8, 5, 5, 5, 4, 4, 2, 2, 2, 0, -1, -3, -4, -4, -7, -7, -8, 10,
                   15, 62, 62, 78, 117, 161, 190, 227, 188, 205, 242, 266, 432, 409, 372, 333, 326, 289, 250, 230,
                   214, 185, 162, 94], float)
REF_LH_T = pd.date_range("1998-06-20 00:30", periods=48, freq="30min")


def load_run(cfg):
    files = sorted(glob.glob(f"{N}/runs/bondville_{cfg}/out/*.LDASOUT_DOMAIN1"))
    parts = []
    for f in files:                                   # plain open_dataset per daily file (no dask needed)
        with xr.open_dataset(f, decode_times=False) as ds:
            t = pd.to_datetime([s.decode() if isinstance(s, bytes) else str(s) for s in ds["Times"].values],
                               format="%Y-%m-%d_%H:%M:%S")
            d = {v: ds[v].values.reshape(len(t), -1)[:, 0] for v in ("LH", "HFX", "TRAD", "TG", "TV", "FVEG", "LAI") if v in ds}
            d["TSLB1"] = ds["SOIL_T"].values[:, 0, 0, 0] if "SOIL_T" in ds else np.full(len(t), np.nan)
            parts.append(pd.DataFrame(d, index=t))
    out = pd.concat(parts).sort_index()
    out = out[~out.index.duplicated()]
    nfill = int((out <= -9000).to_numpy().sum())
    out = out.mask(out <= -9000)                      # HRLDAS writes -9999 fill values (e.g. TRAD at the first output step)
    print(f"[{cfg}] fill values (-9999) masked: {nfill}")
    log = open(f"{N}/runs/bondville_{cfg}/hrldas.log").read()
    rows = re.findall(r"\*\*\*DATE=(\d{4}-\d\d-\d\d_\d\d:\d\d:\d\d)\s+([-\d.]+)\s+([-\d.]+)", log)
    lg = pd.DataFrame({"TSLB1": [float(a) for _, a, _ in rows], "LAI": [float(b) for _, _, b in rows]},
                      index=pd.to_datetime([d for d, _, _ in rows], format="%Y-%m-%d_%H:%M:%S"))
    return out, lg


res = {}
for cfg in CONFIGS:
    out, lg = load_run(cfg)
    res[cfg] = (out, lg)
    d1 = lg["TSLB1"].reindex(REF_T).to_numpy() - REF_TSLB
    lh = out["LH"].reindex(REF_LH_T).to_numpy()
    ok = np.isfinite(lh)
    r = np.corrcoef(lh[ok], REF_LH[ok])[0, 1]
    print(f"\n=== CONFIG={cfg}: {len(out)} half-hours, {out.index[0]} .. {out.index[-1]}")
    print(f"  run-log TSLB1, 02 Jan 00:00-07:30 (16 steps): max |ours - tutorial| = {np.nanmax(np.abs(d1)):.5f} K, "
          f"mean = {np.nanmean(d1):+.5f} K")
    print("    ours    : " + " ".join(f"{v:.3f}" for v in lg['TSLB1'].reindex(REF_T)))
    print("    tutorial: " + " ".join(f"{v:.3f}" for v in REF_TSLB))
    print(f"  LH 20 Jun (48 steps vs digitised ncview): r = {r:.3f}, RMSE = {np.sqrt(np.mean((lh[ok]-REF_LH[ok])**2)):.1f}, "
          f"bias = {np.mean(lh[ok]-REF_LH[ok]):+.1f} W m-2; peak ours {np.nanmax(lh):.0f} at {REF_LH_T[np.nanargmax(lh)]:%H:%M}, "
          f"tutorial {REF_LH.max():.0f} at {REF_LH_T[REF_LH.argmax()]:%H:%M} UTC")
    yr = out.loc["1998-01-02":"1998-12-31"]
    print(f"  year: mean LH {yr['LH'].mean():.1f}, mean HFX {yr['HFX'].mean():.1f} W m-2, ET {yr['LH'].mean()/2.45e6*86400*365:.0f} mm; "
          f"Jun-Aug mean LH {yr.loc['1998-06':'1998-08','LH'].mean():.1f}; TRAD range {yr['TRAD'].min():.1f}-{yr['TRAD'].max():.1f} K")

# Minimal figure (project plot style): tutorial reference (markers) vs our run with the tutorial's options (line).
import matplotlib.dates as mdates, os
ours = res["tutorial"]
fig, ax = plt.subplots(1, 2, figsize=(11, 3.6), layout="constrained")
a = ax[0]
a.plot(REF_LH_T, REF_LH, "o", color="k", ms=4, mfc="none", label="NCAR tutorial")
a.plot(REF_LH_T, ours[0]["LH"].reindex(REF_LH_T), "-", color="tab:red", lw=1.4, label="this run")
a.set_ylabel("LH (W m$^{-2}$)"); a.set_xlabel("20 Jun 1998 (UTC)")
a.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 6))); a.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
a.legend(frameon=False)
a = ax[1]
a.plot(REF_T, REF_TSLB, "o", color="k", ms=4, mfc="none", label="NCAR tutorial")
a.plot(REF_T, ours[1]["TSLB1"].reindex(REF_T), "-", color="tab:red", lw=1.4, label="this run")
a.set_ylabel("Top-layer soil temperature (K)"); a.set_xlabel("2 Jan 1998 (UTC)")
a.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 8, 2))); a.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
for a in ax: a.spines[["top", "right"]].set_visible(False)
os.makedirs(f"{N}/runs/figs", exist_ok=True)
fig.savefig(f"{N}/runs/figs/bondville_reproduce.png", dpi=130); print("\nsaved runs/figs/bondville_reproduce.png")
