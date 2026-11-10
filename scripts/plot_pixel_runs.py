"""Results of the patched Noah-MP forward runs (runs/pixel_<NAME>_2018_2019; 2018 = spin-up), 2019, no MODIS.
Rows: BV, VR, KG. Columns: TRAD (half-hourly), prescribed L4 root-zone theta with the run's s* and s_w, BTRAN (daily mean of
daytime steps, SWFORC > 200 W m-2). Minimal style. Out: runs/figs/pixel_runs_2019.png; prints a short summary per pixel."""
import glob, re
import numpy as np, pandas as pd, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.dates as mdates

N = "/storage/project/r-rbras6-0/jjason3/proposal_drought/noah"; YEAR = 2019; NAMES = ["BV", "VR", "KG"]


def run_series(name):
    rows = []
    for f in sorted(glob.glob(f"{N}/runs/pixel_{name}_2018_2019/out/*.LDASOUT_DOMAIN1")):
        with xr.open_dataset(f, decode_times=False) as ds:
            t = pd.to_datetime([s.decode() if isinstance(s, bytes) else str(s) for s in ds["Times"].values], format="%Y-%m-%d_%H:%M:%S")
            d = {v: ds[v].values.reshape(len(t), -1)[:, 0] for v in ("TRAD", "TV", "TGV", "TGB", "BTRAN", "SWFORC", "LH", "HFX",
                                                                      "ETRAN", "EDIR", "LAI", "FVEG")}
            for k in range(4): d[f"SM{k}"] = ds["SOIL_M"].values[:, 0, k, 0]
            rows.append(pd.DataFrame(d, index=t))
    s = pd.concat(rows).sort_index(); s = s[~s.index.duplicated()]
    return s.mask(s <= -9000)


fig, ax = plt.subplots(len(NAMES), 3, figsize=(15, 3.2 * len(NAMES)), layout="constrained")
for p, name in enumerate(NAMES):
    info = dict(l.strip().split("=", 1) for l in open(f"{N}/runs/pixel_{name}_2018_2019/forcing/pixel_info.txt"))
    sstar, swilt = float(info["LAIO_SSTAR"]), float(info["LAIO_SWILT"])
    M = run_series(name); M = M[M.index.year == YEAR]
    day = M[M.SWFORC > 200]; bd = day.BTRAN.resample("1D").mean()
    theta = (0.1 * M.SM0 + 0.3 * M.SM1 + 0.6 * M.SM2)                 # root zone 0-1 m as the model saw it (layers 2-3 = L4 root zone)
    th_d = M.SM1.resample("1D").mean()                                 # layer 2 = prescribed L4 root zone
    print(f"{name} (IGBP {info['IGBP']}): s* {sstar:.3f}, s_w {swilt:.3f}; L4 root zone {th_d.min():.3f}-{th_d.max():.3f}; "
          f"days below s* {(th_d < sstar).mean():.0%}; daytime BTRAN mean {day.BTRAN.mean():.2f}, min daily {bd.min():.2f}; "
          f"TRAD {M.TRAD.min() - 273.15:.1f} to {M.TRAD.max() - 273.15:.1f} C; mean LAI {M.LAI.mean():.2f}, FVEG {M.FVEG.mean():.2f}; "
          f"LH {M.LH.mean():.0f}, HFX {M.HFX.mean():.0f} W m-2")
    a = ax[p, 0]; a.plot(M.index, M.TRAD - 273.15, "-", color="k", lw=0.3); a.set_ylabel(f"{name}  TRAD ($^\\circ$C)")
    b = ax[p, 1]; b.plot(th_d.index, th_d.values, "-", color="tab:blue", lw=1.2, label="L4 root zone")
    b.axhline(sstar, color="tab:red", lw=1.0, label="s*"); b.axhline(swilt, color="tab:brown", lw=1.0, ls="--", label="$s_w$")
    b.set_ylabel(r"$\theta$ (m$^3$ m$^{-3}$)")
    if p == 0: b.legend(frameon=False, fontsize=8)
    c = ax[p, 2]; c.plot(bd.index, bd.values, "-", color="tab:green", lw=1.2); c.set_ylabel("BTRAN"); c.set_ylim(-0.02, 1.05)
for a in ax.ravel():
    a.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 4, 7, 10])); a.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    a.spines[["top", "right"]].set_visible(False)
fig.savefig(f"{N}/runs/figs/pixel_runs_{YEAR}.png", dpi=120); print(f"saved runs/figs/pixel_runs_{YEAR}.png")
