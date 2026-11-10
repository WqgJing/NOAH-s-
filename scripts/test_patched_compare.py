"""Checks for the patched Noah-MP (hpc/test_patched_bondville.slurm), Bondville 1998.
  1  R  == stock run, every output variable, bit for bit (switches off must reproduce NCAR).
  2  B4 : BTRAN(t) == Laio(theta_root) with theta_root = thickness-weighted mean of SOIL_W over the 4 root layers at the start
          of the step (= output of the previous step); exact up to float32 output rounding.
  3  P  : BTRAN(t) == Laio(SMROOT) of the prescribed forcing (best of lags -1/0/+1 step); SOIL_M layer 1 follows SMSURF.
  4  S1/S2/S3 vs R: direction of the nuisance scales in Jun-Aug daytime (SWFORC > 200 W m-2).
  5  physics time (SFLXtot) vs total run time."""
import glob, os, re
import numpy as np, pandas as pd, xarray as xr

N = "/storage/project/r-rbras6-0/jjason3/proposal_drought/noah"; T = f"{N}/runs/test_patched"
W = np.array([0.1, 0.3, 0.6]) / 1.0                          # root layers of cropland, MODIS table NROOT(12) = 3: thickness / 1 m
SSTAR, SWILT = 0.387, 0.120
laio = lambda th: np.maximum(1e-6, np.clip((th - SWILT) / (SSTAR - SWILT), 0.0, 1.0))


def times(ds):
    return pd.to_datetime([s.decode() if isinstance(s, bytes) else str(s) for s in ds["Times"].values], format="%Y-%m-%d_%H:%M:%S")


def series(d):
    rows = []
    for f in sorted(glob.glob(f"{d}/*.LDASOUT_DOMAIN1")):
        with xr.open_dataset(f, decode_times=False) as ds:
            t = times(ds); dd = {v: ds[v].values.reshape(len(t), -1)[:, 0] for v in
                                 ("LH", "HFX", "TRAD", "EDIR", "ETRAN", "SWFORC", "BTRAN", "T2MV") if v in ds}
            for k in range(4):
                dd[f"SW{k}"] = ds["SOIL_W"].values[:, 0, k, 0]; dd[f"SM{k}"] = ds["SOIL_M"].values[:, 0, k, 0]
            rows.append(pd.DataFrame(dd, index=t))
    s = pd.concat(rows).sort_index(); return s[~s.index.duplicated()]


# 1 regression: R vs stock, bit for bit
stock_files = sorted(glob.glob(f"{N}/runs/bondville_shipped/out/*.LDASOUT_DOMAIN1"))
nvar = nbad = 0; bad = set()
for fs in stock_files:
    fr = f"{T}/R/out/{os.path.basename(fs)}"
    with xr.open_dataset(fs, decode_times=False) as a, xr.open_dataset(fr, decode_times=False) as b:
        for v in a.data_vars:
            nvar += 1
            if not np.array_equal(a[v].values, b[v].values, equal_nan=(a[v].dtype.kind == "f")):
                nbad += 1; bad.add(v)
        has_btran = "BTRAN" in b
print(f"[1] R vs stock: {len(stock_files)} daily files, {nvar} variable-files compared, mismatches {nbad} {sorted(bad)} "
      f"-> {'PASS (bit-identical)' if nbad == 0 else 'FAIL'}; BTRAN in patched output: {has_btran}")

R, B4, P = series(f"{T}/R/out"), series(f"{T}/B4/out"), series(f"{T}/P/out")

# 2 Laio option 4 against its formula on the start-of-step state
th0 = (B4[[f"SW{k}" for k in range(3)]].to_numpy() * W).sum(1)   # end-of-step state at t (3 root layers)
pred = pd.Series(laio(th0), index=B4.index).shift(1)              # -> start-of-step state of the next step
ok = B4["BTRAN"].notna() & (B4["BTRAN"] >= 0) & pred.notna()
d2 = (B4["BTRAN"] - pred)[ok].abs()
print(f"[2] B4 BTRAN vs Laio(start-of-step root-zone liquid water): n {ok.sum()}, max |diff| {d2.max():.2e}, "
      f"99th pct {d2.quantile(0.99):.2e}, share < 1e-4: {(d2 < 1e-4).mean():.1%} -> {'PASS' if d2.quantile(0.99) < 1e-4 else 'CHECK'}; "
      f"BTRAN range {B4.BTRAN[ok].min():.3f}-{B4.BTRAN[ok].max():.3f}, mean R(opt 1) {R.BTRAN[R.BTRAN >= 0].mean():.3f} vs B4 {B4.BTRAN[ok].mean():.3f}")

# 3 prescription
S = series(f"{N}/runs/bondville_shipped/out")
smroot = (0.1 * S.SM0 + 0.3 * S.SM1 + 0.6 * S.SM2) / 1.0; smsurf = S.SM0
best = None
for lag in (-1, 0, 1):
    pr = pd.Series(laio(smroot.to_numpy()), index=smroot.index).shift(lag).reindex(P.index)
    okp = P["BTRAN"].notna() & (P["BTRAN"] >= 0) & pr.notna(); dd = (P["BTRAN"] - pr)[okp].abs()
    print(f"    P lag {lag:+d}: BTRAN vs Laio(SMROOT) max |diff| {dd.max():.2e}, median {dd.median():.2e}")
    if best is None or dd.median() < best[1]: best = (lag, dd.median(), dd.max(), dd.quantile(0.99))
lag = best[0]
s1 = (P.SM0 - smsurf.shift(lag).reindex(P.index)).abs()
r2 = (P.SM1 - smroot.shift(lag).reindex(P.index)).abs()
print(f"[3] P: best lag {lag:+d} step; BTRAN vs Laio(SMROOT) 99th pct |diff| {best[3]:.2e} -> {'PASS' if best[3] < 1e-4 else 'CHECK'}; "
      f"end-of-step SOIL_M layer 1 vs SMSURF: median |diff| {s1.median():.4f}, 99th {s1.quantile(0.99):.4f}; "
      f"layer 2 vs SMROOT: median {r2.median():.4f}, 99th {r2.quantile(0.99):.4f}")

# 4 nuisance directions (Jun-Aug daytime)
def summer_day(df):
    m = (df.index.month >= 6) & (df.index.month <= 8) & (df.SWFORC > 200); return df[m]
r = summer_day(R)
for case, var, sign, txt in (("S1", "ETRAN", +1, "conductance x2 -> more transpiration"),
                             ("S1", "TRAD", -1, "conductance x2 -> cooler surface"),
                             ("S2", "EDIR", -1, "soil resistance x2 -> less soil evaporation"),
                             ("S3", "TRAD", -1, "roughness x2 -> daytime surface closer to air (cooler)")):
    s = summer_day(series(f"{T}/{case}/out")); dv = s[var].mean() - r[var].mean()
    unit = "mm/s" if var in ("ETRAN", "EDIR") else "K"
    print(f"[4] {case} {txt}: mean {var} {r[var].mean():.4g} -> {s[var].mean():.4g} ({dv:+.3g} {unit}) "
          f"-> {'PASS' if np.sign(dv) == sign else 'CHECK'}")

# 5 timing
for c in ("R", "B4", "P", "S1", "S2", "S3"):
    m = re.findall(r"Cumulative:\s+([\d.]+)\s+SFLX:\s+[\d.]+\s+SFLXtot:\s+([\d.]+)", open(f"{T}/{c}/hrldas.log").read())
    if m: print(f"[5] {c}: total {float(m[-1][0]):.1f} s, physics (SFLX) {float(m[-1][1]):.1f} s per point-year")
