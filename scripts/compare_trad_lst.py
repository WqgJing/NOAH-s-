"""Noah-MP TRAD (patched model, our inputs; runs/pixel_<NAME>_2018_2019, 2018 = spin-up) vs MODIS LST at the test pixels, 2019.
Matching: a MODIS day/night value is used when >= 50 % of the 0.25 deg cell's 1 km pixels are clear (QC err2); its time is the mean
view time (local solar hour) converted to UTC with the pixel longitude; TRAD (half-hourly, end-of-step) is interpolated linearly to it.
Prints n, bias (TRAD - LST), RMSE and r per pixel, satellite and pass; also BTRAN use (share of daytime steps with BTRAN < 1) and the
canopy / ground temperatures at the Aqua day pass. Figure (minimal style): runs/figs/trad_vs_lst_2019.png."""
import glob
import numpy as np, pandas as pd, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.dates as mdates

N = "/storage/project/r-rbras6-0/jjason3/proposal_drought/noah"; YEAR = 2019; MINCLEAR = 0.5
S = np.load(f"{N}/data/lst_points.npz", allow_pickle=True); L = np.load(f"{N}/data/l4_points.npz", allow_pickle=True)
dates = pd.DatetimeIndex(S["date"]); sats = list(S["sat"]); names = list(S["name"])


def run_series(name):
    rows = []
    for f in sorted(glob.glob(f"{N}/runs/pixel_{name}_2018_2019/out/*.LDASOUT_DOMAIN1")):
        with xr.open_dataset(f, decode_times=False) as ds:
            t = pd.to_datetime([s.decode() if isinstance(s, bytes) else str(s) for s in ds["Times"].values], format="%Y-%m-%d_%H:%M:%S")
            rows.append(pd.DataFrame({v: ds[v].values.reshape(len(t), -1)[:, 0] for v in
                                      ("TRAD", "TV", "TGV", "TGB", "FVEG", "LAI", "BTRAN", "SWFORC", "T2MV", "LH")}, index=t))
    s = pd.concat(rows).sort_index(); s = s[~s.index.duplicated()]
    return s.mask(s <= -9000)


fig, ax = plt.subplots(len(names), 2, figsize=(13, 3.4 * len(names)), layout="constrained", gridspec_kw={"width_ratios": [3, 1]})
for p, name in enumerate(names):
    M = run_series(name); lon = float(L["lon"][list(L["name"]).index(name)])
    tsec = M.index.values.astype("datetime64[s]").astype(np.int64)
    day = (M.index.year == YEAR) & (M.SWFORC > 200)
    print(f"\n== {name} ({YEAR}): mean LAI {M.LAI[M.index.year == YEAR].mean():.2f}, FVEG {M.FVEG[M.index.year == YEAR].mean():.2f}; "
          f"daytime steps with BTRAN < 1: {(M.BTRAN[day] < 0.999).mean():.0%}, BTRAN min {M.BTRAN[day].min():.2f}")
    keep = {}
    for s_i, sat in enumerate(sats):
        for pas, kl, kf, kv in (("day", "ld", "fd", "vtd"), ("night", "ln", "fn", "vtn")):
            lst, fr, vt = S[kl][s_i, :, p], S[kf][s_i, :, p], S[kv][s_i, :, p]
            ok = (dates.year == YEAR) & np.isfinite(lst) & (fr >= MINCLEAR) & np.isfinite(vt)
            tu = dates[ok] + pd.to_timedelta(vt[ok] - lon / 15.0, "h")
            trad = np.interp(tu.values.astype("datetime64[s]").astype(np.int64), tsec, M.TRAD.to_numpy())
            d = trad - lst[ok]
            r = np.corrcoef(trad, lst[ok])[0, 1] if ok.sum() > 2 else np.nan
            print(f"  {sat:5s} {pas:5s}: n {ok.sum():3d}, bias {np.mean(d):+.2f} K, RMSE {np.sqrt(np.mean(d ** 2)):.2f} K, r {r:.3f}")
            keep[(sat, pas)] = (tu, lst[ok], trad)
    tu, lst, trad = keep[("aqua", "day")]
    tv = np.interp(tu.values.astype("datetime64[s]").astype(np.int64), tsec, M.TV.to_numpy())
    tgb = np.interp(tu.values.astype("datetime64[s]").astype(np.int64), tsec, M.TGB.to_numpy())
    print(f"  at Aqua day pass: mean TRAD {np.mean(trad) - 273.15:.1f} C, canopy TV {np.mean(tv) - 273.15:.1f} C, "
          f"bare ground TGB {np.mean(tgb) - 273.15:.1f} C, MODIS {np.mean(lst) - 273.15:.1f} C")
    a = ax[p, 0]
    a.plot(tu, lst - 273.15, "o", ms=3, color="k", mfc="none", label="MODIS Aqua day LST")
    a.plot(tu, trad - 273.15, "-", lw=1.0, color="tab:red", label="Noah-MP TRAD")
    a.set_ylabel(f"{name}  temperature ($^\\circ$C)")
    a.xaxis.set_major_locator(mdates.MonthLocator()); a.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    if p == 0: a.legend(frameon=False)
    b = ax[p, 1]
    for (sat, pas), c in ((("aqua", "day"), "tab:red"), (("terra", "day"), "tab:orange"), (("aqua", "night"), "tab:blue"), (("terra", "night"), "tab:cyan")):
        tu_, l_, t_ = keep[(sat, pas)]; b.plot(l_ - 273.15, t_ - 273.15, "o", ms=2.5, color=c, label=f"{sat} {pas}")
    lim = [min(b.get_xlim()[0], b.get_ylim()[0]), max(b.get_xlim()[1], b.get_ylim()[1])]
    b.plot(lim, lim, "-", color="0.5", lw=0.8); b.set_xlim(lim); b.set_ylim(lim)
    b.set_xlabel("MODIS LST ($^\\circ$C)"); b.set_ylabel("TRAD ($^\\circ$C)")
    if p == 0: b.legend(frameon=False, fontsize=8)
for a in ax.ravel(): a.spines[["top", "right"]].set_visible(False)
fig.savefig(f"{N}/runs/figs/trad_vs_lst_{YEAR}.png", dpi=120); print(f"\nsaved runs/figs/trad_vs_lst_{YEAR}.png")
