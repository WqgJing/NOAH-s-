"""Python port (pyport) vs the Fortran Noah-MP reference runs, 1 Jun - 30 Sep 2018 (runs/pyref_*: VR, KG with reference s*;
VRs, KGsn with s* = 0.25 so the Laio stress is active). Python output: runs/<run>/pyport_out.npz (scripts/pyport_compare.py);
Fortran output: runs/<run>/out/*.LDASOUT_DOMAIN1. Rows = runs; columns: TRAD (one week), LH (same week), BTRAN (full window),
TRAD difference Python - Fortran (full window). Minimal style. Out: runs/figs/pyport_vs_fortran.png"""
import glob
import numpy as np, pandas as pd, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.dates as mdates

N = "/storage/project/r-rbras6-0/jjason3/proposal_drought/noah"
RUNS = ["pyref_VR", "pyref_KG", "pyref_VRs", "pyref_KGsn"]
LABEL = {"pyref_VR": "VR", "pyref_KG": "KG", "pyref_VRs": "VR, s* 0.25", "pyref_KGsn": "KG, s* 0.25"}
WEEK = (pd.Timestamp("2018-07-01"), pd.Timestamp("2018-07-08"))
VARS = ("TRAD", "LH", "BTRAN")


def fortran(run):
    parts = []
    for f in sorted(glob.glob(f"{N}/runs/{run}/out/*.LDASOUT_DOMAIN1")):
        with xr.open_dataset(f, decode_times=False) as ds:
            t = pd.to_datetime([s.decode() if isinstance(s, bytes) else str(s) for s in ds["Times"].values], format="%Y-%m-%d_%H:%M:%S")
            parts.append(pd.DataFrame({v: ds[v].values.reshape(len(t), -1)[:, 0] for v in VARS}, index=t))
    s = pd.concat(parts).sort_index(); return s[~s.index.duplicated()].mask(lambda d: d <= -9000)


fig, ax = plt.subplots(len(RUNS), 4, figsize=(17, 3.0 * len(RUNS)), layout="constrained")
for r, run in enumerate(RUNS):
    Z = np.load(f"{N}/runs/{run}/pyport_out.npz")
    P = pd.DataFrame({v: Z[v] for v in VARS}, index=pd.DatetimeIndex(Z["time"]))
    F = fortran(run).reindex(P.index)
    w = (P.index >= WEEK[0]) & (P.index < WEEK[1])
    d = P.TRAD - F.TRAD
    print(f"{run}: n {len(P)}; TRAD diff max {d.abs().max():.2e} K, RMS {np.sqrt((d ** 2).mean()):.2e} K; "
          f"LH diff max {(P.LH - F.LH).abs().max():.2e}; BTRAN diff max {(P.BTRAN - F.BTRAN).abs().max():.2e}")
    a = ax[r, 0]
    a.plot(F.index[w], F.TRAD[w] - 273.15, "-", color="k", lw=2.2, label="Fortran")
    a.plot(P.index[w], P.TRAD[w] - 273.15, "--", color="tab:red", lw=1.2, label="Python")
    a.set_ylabel(f"{LABEL[run]}\nTRAD ($^\\circ$C)")
    a = ax[r, 1]
    a.plot(F.index[w], F.LH[w], "-", color="k", lw=2.2); a.plot(P.index[w], P.LH[w], "--", color="tab:red", lw=1.2)
    a.set_ylabel("LH (W m$^{-2}$)")
    for a in ax[r, :2]:
        a.xaxis.set_major_locator(mdates.DayLocator(interval=2)); a.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    a = ax[r, 2]
    a.plot(F.index, F.BTRAN, "-", color="k", lw=2.2); a.plot(P.index, P.BTRAN, "--", color="tab:red", lw=1.2)
    a.set_ylabel("BTRAN"); a.set_ylim(-0.02, 1.05)
    a = ax[r, 3]
    a.plot(d.index, 1000 * d, "-", color="tab:blue", lw=0.5); a.axhline(0, color="0.6", lw=0.8)
    a.set_ylabel("TRAD Python $-$ Fortran (mK)")
    for a in ax[r, 2:]:
        a.xaxis.set_major_locator(mdates.MonthLocator()); a.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
ax[0, 0].legend(frameon=False)
for a in ax.ravel(): a.spines[["top", "right"]].set_visible(False)
fig.savefig(f"{N}/runs/figs/pyport_vs_fortran.png", dpi=110); print("saved runs/figs/pyport_vs_fortran.png")
