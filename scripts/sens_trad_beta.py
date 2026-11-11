"""Sensitivity of simulated TRAD to the water-stress factor beta (BTRAN), Python port, our full configuration,
VR and KG, 1 Jun - 30 Sep 2018 (runs/pyref_VR, runs/pyref_KG: same inputs as the verified Fortran reference runs).
  E1  fixed beta in {1, 0.8, 0.6, 0.4, 0.2, 0.05}, ignoring soil moisture (pyport hook p.BtranOverride): pure closure response.
  E2  s* scan: s* = quantile q of the prescribed L4 root-zone theta in the window, q in {0, 0.25, 0.5, 0.75, 1}
      (root zone below s* for a share q of the time), s_w fixed at the run's value; beta from the Laio function as in the model.
Metric: TRAD at the Aqua day overpass (13:30 local solar = 13.5 - lon/15 h UTC), linear in time between half-hourly outputs;
also canopy (TV) and bare-ground (TGB) temperature and transpiration at that time. Differences are vs beta = 1.
Check: the default run (no override) must equal the verified port output runs/<run>/pyport_out.npz bit for bit.
Out: printed tables, runs/figs/sens_trad_beta.png (minimal style), runs/sens_trad_beta.npz"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.dates as mdates

N = "/storage/project/r-rbras6-0/jjason3/proposal_drought/noah"; sys.path.insert(0, N)
from pyport import PointModel

RUNS = {"VR": f"{N}/runs/pyref_VR", "KG": f"{N}/runs/pyref_KG"}
BETAS = [1.0, 0.8, 0.6, 0.4, 0.2, 0.05]
QS = [0.0, 0.25, 0.5, 0.75, 1.0]
VARS = ("TRAD", "TV", "TGB", "ETRAN", "LH", "BTRAN")


def run(rd, beta=None, sstar=None):
    m = PointModel(rd)
    if beta is not None: m.p.BtranOverride = beta
    if sstar is not None: m.p.SoilMoistureLaioStar = float(sstar)
    t, o = m.run()
    return pd.DataFrame({v: o[v] for v in VARS}, index=pd.DatetimeIndex(t)), m


def at_overpass(df, lon):
    days = pd.date_range(df.index[0].normalize() + pd.Timedelta("1D"), df.index[-1].normalize() - pd.Timedelta("1D"), freq="D")
    tu = days + pd.to_timedelta(13.5 - lon / 15.0, "h")
    x = df.index.values.astype("datetime64[s]").astype(np.int64); xi = tu.values.astype("datetime64[s]").astype(np.int64)
    return pd.DataFrame({v: np.interp(xi, x, df[v].to_numpy()) for v in df}, index=days)


RES = {}
fig, ax = plt.subplots(len(RUNS), 3, figsize=(16, 3.6 * len(RUNS)), layout="constrained")
for r, (name, rd) in enumerate(RUNS.items()):
    base, m0 = run(rd)
    ref = np.load(f"{rd}/pyport_out.npz")
    print(f"\n== {name}: default run vs verified port output, TRAD max |diff| {np.max(np.abs(base.TRAD.to_numpy() - ref['TRAD'])):.1e} K "
          f"(must be 0); BTRAN range {base.BTRAN.min():.3f}-{base.BTRAN.max():.3f}; s* {m0.p.SoilMoistureLaioStar:.3f}, s_w {m0.p.SoilMoistureLaioWilt:.3f}")
    F = np.load(f"{rd}/forcing/forcing_arrays.npz"); lon = float(F["lon"])
    ft = pd.DatetimeIndex(F["time"]); win = (ft >= base.index[0]) & (ft <= base.index[-1])
    th = F["SMROOT"][win]
    # E1 fixed beta
    e1 = {b: at_overpass(run(rd, beta=b)[0], lon) for b in BETAS}
    b1 = e1[1.0]
    print(f"E1 fixed beta (13:30 local, Jun-Sep mean; differences vs beta = 1):")
    print(f"  {'beta':>5} {'TRAD':>8} {'dTRAD':>7} {'dTV':>7} {'dTGB':>7} {'ETRAN mm/d':>11}")
    rows = []
    for b, d in e1.items():
        rows.append((b, (d.TRAD - b1.TRAD).mean(), (d.TV - b1.TV).mean(), (d.TGB - b1.TGB).mean()))
        print(f"  {b:5.2f} {d.TRAD.mean() - 273.15:8.2f} {rows[-1][1]:+7.2f} {rows[-1][2]:+7.2f} {rows[-1][3]:+7.2f} {d.ETRAN.mean() * 86400:11.3f}")
    rows = np.array(rows)
    # E2 s* scan
    swilt = m0.p.SoilMoistureLaioWilt
    e2 = {}
    print(f"E2 s* scan (root-zone theta in window {th.min():.3f}-{th.max():.3f}):")
    for q in QS:
        s = max(float(np.quantile(th, q)), swilt + 0.01)
        dfq, _ = run(rd, sstar=s); d = at_overpass(dfq, lon); e2[q] = (s, d)
        print(f"  q {q:4.2f}: s* {s:.3f}, daytime-overpass BTRAN mean {d.BTRAN.mean():.2f} (min {d.BTRAN.min():.2f}), "
              f"dTRAD vs beta=1 mean {(d.TRAD - b1.TRAD).mean():+.2f} K, max {(d.TRAD - b1.TRAD).max():+.2f} K")
    RES[name] = dict(e1=rows, e2=[(q, s, (d.TRAD - b1.TRAD).to_numpy()) for q, (s, d) in e2.items()])
    a = ax[r, 0]
    a.plot(rows[:, 0], rows[:, 1], "o-", color="k", lw=1.4, label="TRAD")
    a.plot(rows[:, 0], rows[:, 2], "s--", color="tab:green", lw=1.2, label="canopy TV")
    a.plot(rows[:, 0], rows[:, 3], "^:", color="tab:brown", lw=1.2, label="bare ground TGB")
    a.set_xlabel(r"$\beta$ (fixed)"); a.set_ylabel(f"{name}\n$\\Delta$T at 13:30 vs $\\beta$=1 (K)"); a.invert_xaxis()
    if r == 0: a.legend(frameon=False)
    a = ax[r, 1]
    for b, c in zip(BETAS, plt.cm.viridis(np.linspace(0, 0.9, len(BETAS)))):
        a.plot(e1[b].index, e1[b].TRAD - 273.15, "-", color=c, lw=1.0, label=f"$\\beta$ = {b}")
    a.set_ylabel(r"TRAD at 13:30 ($^\circ$C)")
    if r == 0: a.legend(frameon=False, fontsize=8, ncol=2)
    a = ax[r, 2]
    for (q, (s, d)), c in zip(e2.items(), plt.cm.plasma(np.linspace(0, 0.85, len(QS)))):
        a.plot(d.index, d.TRAD - b1.TRAD, "-", color=c, lw=1.0, label=f"s* {s:.3f} (below {q:.0%})")
    a.axhline(0, color="0.6", lw=0.8); a.set_ylabel(r"$\Delta$TRAD at 13:30 vs $\beta$=1 (K)")
    a.legend(frameon=False, fontsize=8)
    for a in ax[r, 1:]:
        a.xaxis.set_major_locator(mdates.MonthLocator()); a.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
for a in ax.ravel(): a.spines[["top", "right"]].set_visible(False)
fig.savefig(f"{N}/runs/figs/sens_trad_beta.png", dpi=110)
np.savez(f"{N}/runs/sens_trad_beta.npz", **{f"{k}_e1": v["e1"] for k, v in RES.items()})
print("\nsaved runs/figs/sens_trad_beta.png")
