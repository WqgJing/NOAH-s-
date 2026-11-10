"""Bondville example (shipped options): half-hourly TRAD for 1998. Minimal figure -> runs/figs/bondville_trad.png"""
import glob
import pandas as pd, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.dates as mdates

N = "/storage/project/r-rbras6-0/jjason3/proposal_drought/noah"
parts = []
for f in sorted(glob.glob(f"{N}/runs/bondville_shipped/out/*.LDASOUT_DOMAIN1")):
    with xr.open_dataset(f, decode_times=False) as ds:
        t = pd.to_datetime([s.decode() if isinstance(s, bytes) else str(s) for s in ds["Times"].values], format="%Y-%m-%d_%H:%M:%S")
        parts.append(pd.Series(ds["TRAD"].values.reshape(len(t), -1)[:, 0], index=t))
trad = pd.concat(parts).sort_index()
trad = trad[~trad.index.duplicated()].mask(lambda s: s <= -9000)

fig, ax = plt.subplots(figsize=(10, 3.2), layout="constrained")
ax.plot(trad.index, trad.values, "-", color="k", lw=0.4)
ax.set_ylabel("TRAD (K)")
ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
ax.set_xlim(trad.index[0], trad.index[-1])
ax.spines[["top", "right"]].set_visible(False)
fig.savefig(f"{N}/runs/figs/bondville_trad.png", dpi=130)
print("saved runs/figs/bondville_trad.png", len(trad), "values")
