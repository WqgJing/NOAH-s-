"""Compare the LOCAL SMAP L4 daily root zone (sbi_train/runs/l4/l4_rootzone_<year>.npy, the product the Laio map uses) with the
GEE SPL4SMGP/008 3-hourly root zone fetched for the Noah-MP test pixels (noah/data/l4_points.npz), as daily means (UTC days).
Prints per pixel and year: n days, mean local, mean difference (GEE - local), max |difference|, r."""
import numpy as np, pandas as pd

D = "/storage/project/r-rbras6-0/jjason3/proposal_drought"
L = np.load(f"{D}/noah/data/l4_points.npz", allow_pickle=True)
B = np.load(f"{D}/sbi_train/runs/vegsbi/bank.npz")
grid = pd.read_csv("/storage/project/r-rbras6-0/jjason3/forward_sensitivity_lib/data/global/global_grid_025deg.csv").sort_values("grid_id")
row_of = {g: i for i, g in enumerate(grid.grid_id.to_numpy().astype(np.int64))}
t = pd.DatetimeIndex(L["time"])
for j, name in enumerate(L["name"]):
    b = int(L["bank"][j]); gid = int(B["grid_id"][b]); r = row_of[gid]
    gee = pd.Series(L["sm_rootzone"][:, j], index=t).resample("1D").mean()
    g = grid.iloc[r]
    print(f"{name}: bank {b}, grid_id {gid} at ({g.lat:.3f}, {g.lon:.3f}); pixel ({float(L['lat'][j]):.3f}, {float(L['lon'][j]):.3f})")
    for y in range(2018, 2025):
        loc = np.load(f"{D}/sbi_train/runs/l4/l4_rootzone_{y}.npy", mmap_mode="r")[r]
        ls = pd.Series(np.asarray(loc, float), index=pd.date_range(f"{y}-01-01", periods=len(loc), freq="D"))
        both = pd.concat([ls, gee], axis=1, keys=["loc", "gee"]).dropna()
        d = both.gee - both["loc"]
        print(f"  {y}: n {len(both):3d}, local mean {both['loc'].mean():.4f}, GEE-local mean {d.mean():+.5f}, max |d| {d.abs().max():.5f}, "
              f"r {np.corrcoef(both['loc'], both.gee)[0, 1]:.5f}")
