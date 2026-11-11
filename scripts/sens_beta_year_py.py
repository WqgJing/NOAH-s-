"""Where / when does TRAD respond to beta? Python port (pyport), our full configuration, VR and KG, 2018-2019 inputs of the
Fortran pixel runs (runs/pixel_<NAME>_2018_2019; 2018 = spin-up, 2019 analysed).
1. Check: default port run vs the Fortran output of the same run, monthly RMS of TRAD (covers winter / frost nights, which the
   Jun-Sep reference runs did not). If the port meets something it does not support (snow on the ground), it stops there and
   the analysis uses the period before that.
2. Fixed beta sweep (p.BtranOverride) b in {1, 0.5, 0.2, 0.05}: TRAD at the Aqua day overpass (13:30 local solar), daily 2019;
   dTRAD vs b = 1; transpiration heat TR and LAI at 13:30 for b = 1 (the "available" closure signal).
Window (2026-11-10, user: "check several months like growing season"): 2019-03-01 .. 10-31, cold start, first 10 days dropped.
Out: printed monthly tables, runs/figs/sens_beta_year_py.png (minimal style)."""
import glob, sys
from datetime import datetime
import numpy as np, pandas as pd, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.dates as mdates

N = "/storage/project/r-rbras6-0/jjason3/proposal_drought/noah"; sys.path.insert(0, N)
from pyport import PointModel
from pyport.driver import ldasout_record

NAMES = ["VR", "KG"]; BETAS = [1.0, 0.5, 0.2, 0.05]
KEEP = ("TRAD", "TV", "TGB", "TR", "LAI", "FVEG", "BTRAN", "ETRAN")


START, END, SPIN = datetime(2019, 3, 1), datetime(2019, 11, 1), 10      # growing-season window, cold start, spin-up days dropped


def run(rd, beta=None):
    m = PointModel(rd)
    m.start = START; m.ntime = int((END - START).total_seconds() / m.dt)   # cold start from the setup file at START
    if beta is not None: m.p.BtranOverride = beta
    m.read_forcing(0); times, rows, stop = [], [], None
    for it in range(1, m.ntime + 1):
        try:
            when, v = m.step(it)
        except NotImplementedError as e:
            stop = (str(times[-1]) if times else "start", str(e)[:120]); break
        times.append(when); r = ldasout_record(v, m.io); rows.append([r[k] for k in KEEP])
    df = pd.DataFrame(np.array(rows, float), columns=KEEP, index=pd.DatetimeIndex(np.array(times, dtype="datetime64[s]")))
    return df.mask(df <= -9000), stop


def fortran(rd):
    parts = []
    for f in sorted(glob.glob(f"{rd}/out/*.LDASOUT_DOMAIN1")):
        with xr.open_dataset(f, decode_times=False) as ds:
            t = pd.to_datetime([s.decode() if isinstance(s, bytes) else str(s) for s in ds["Times"].values], format="%Y-%m-%d_%H:%M:%S")
            parts.append(pd.DataFrame({"TRAD": ds["TRAD"].values.reshape(len(t), -1)[:, 0]}, index=t))
    s = pd.concat(parts).sort_index(); return s[~s.index.duplicated()].mask(lambda d: d <= -9000)


def at_overpass(df, lon):
    days = pd.date_range(START + pd.Timedelta(days=SPIN), END - pd.Timedelta(days=1), freq="D")
    tu = days + pd.to_timedelta(13.5 - lon / 15.0, "h")
    x = df.index.values.astype("datetime64[s]").astype(np.int64); xi = tu.values.astype("datetime64[s]").astype(np.int64)
    ok = (xi >= x[0]) & (xi <= x[-1])
    out = pd.DataFrame({v: np.interp(xi, x, df[v].to_numpy()) for v in df}, index=days)
    return out[ok]


fig, ax = plt.subplots(len(NAMES), 3, figsize=(17, 3.8 * len(NAMES)), layout="constrained")
for r, name in enumerate(NAMES):
    rd = f"{N}/runs/pixel_{name}_2018_2019"; lon = float(np.load(f"{rd}/forcing/forcing_arrays.npz")["lon"])
    sw = {}
    for b in BETAS:
        df, stop = run(rd, beta=b); sw[b] = at_overpass(df, lon)
        print(f"   {name} beta {b}: ran {df.index[0]} .. {df.index[-1]}; stopped early: {stop}")
    b1 = sw[1.0]
    dT = {b: sw[b].TRAD - b1.TRAD for b in BETAS}
    m1 = b1.resample("MS").mean(); mdt = {b: dT[b].resample("MS").mean() for b in BETAS}
    print("   month   LAI  FVEG  TR13:30(W/m2)  dTRAD b=0.5  b=0.2  b=0.05 (K, 13:30)")
    for t in m1.index:
        print(f"   {t:%b}  {m1.LAI[t]:4.2f}  {m1.FVEG[t]:4.2f}  {m1.TR[t]:10.0f}     {mdt[0.5][t]:+6.2f} {mdt[0.2][t]:+6.2f} {mdt[0.05][t]:+6.2f}")
    a = ax[r, 0]
    for b, c in zip(BETAS[1:], ("tab:orange", "tab:red", "darkred")):
        a.plot(dT[b].index, dT[b].values, "-", color=c, lw=1.0, label=f"$\\beta$ = {b}")
    a.axhline(0, color="0.6", lw=0.8); a.set_ylabel(f"{name}\n$\\Delta$TRAD at 13:30 vs $\\beta$=1 (K)")
    if r == 0: a.legend(frameon=False)
    a = ax[r, 1]
    a.plot(b1.index, b1.TR, "-", color="tab:green", lw=1.0); a.set_ylabel(r"transpiration heat at 13:30, $\beta$=1 (W m$^{-2}$)")
    for a in ax[r, :2]:
        a.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 3, 5, 7, 9, 11])); a.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    a = ax[r, 2]
    a.plot(b1.TR, dT[0.05], "o", ms=2.5, color="darkred", mfc="none")
    a.set_xlabel(r"transpiration heat at 13:30, $\beta$=1 (W m$^{-2}$)"); a.set_ylabel(r"$\Delta$TRAD, $\beta$=0.05 (K)")
for a in ax.ravel(): a.spines[["top", "right"]].set_visible(False)
fig.savefig(f"{N}/runs/figs/sens_beta_year_py.png", dpi=110); print("\nsaved runs/figs/sens_beta_year_py.png")
