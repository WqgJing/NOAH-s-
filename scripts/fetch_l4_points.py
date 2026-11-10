"""SMAP L4 (SPL4SMGP, one version for all years) 3-hourly sm_surface + sm_rootzone at a few Noah-MP test pixels, 2018-2024.
Test pixels = the nearest pixel of the CliMA/ERA5 global pixel set (clima/forcing_global/pixels.npz, ERA5 grid points) to:
  BV Bondville IL (cropland; the NCAR example site), VR Vaira Ranch CA (Mediterranean grassland, deep summer drydown),
  KG Kendall AZ (semi-arid grassland).
Out: noah/data/l4_points.npz (time [UTC], sm_surface, sm_rootzone [n_time, n_pix], pixel metadata incl. bank index and
L4 Campbell soil parameters from the bank). Run via sbatch (hpc/fetch_l4_points.slurm)."""
import json, time
import numpy as np, pandas as pd, ee

N = "/storage/project/r-rbras6-0/jjason3/proposal_drought"
KEY = "/storage/project/r-rbras6-0/jjason3/drought/gldasnsmap-dee8cada466e.json"
ASSET = "NASA/SMAP/SPL4SMGP/008"
YEARS = range(2018, 2025)
TARGETS = {"BV": (40.01, -88.37), "VR": (38.41, -120.95), "KG": (31.74, -109.94)}

P = np.load(f"{N}/clima/forcing_global/pixels.npz"); plat, plon, pbank = P["lat"], P["lon"], P["bank_idx"]
B = np.load(f"{N}/sbi_train/runs/vegsbi/bank.npz")
LC = np.load(f"{N}/sbi_train/runs/psi/lc_igbp.npz")["lc"]
sel = {}
for k, (la, lo) in TARGETS.items():
    j = int(np.argmin((plat - la) ** 2 + ((plon - lo) * np.cos(np.radians(la))) ** 2))
    sel[k] = j
    b = int(pbank[j])
    print(f"{k}: target ({la}, {lo}) -> pixel {j} at ({plat[j]:.3f}, {plon[j]:.3f}), bank {b}, IGBP {int(LC[b])}, "
          f"poros {B['poros'][b]:.3f} bee {B['bee'][b]:.2f} psis {B['psis'][b]:.4f}", flush=True)

info = json.load(open(KEY))
ee.Initialize(ee.ServiceAccountCredentials(info["client_email"], KEY), project=info["project_id"],
              opt_url="https://earthengine-highvolume.googleapis.com")
times, SURF, ROOT = None, [], []
for k, j in sel.items():
    pt = ee.Geometry.Point([float(plon[j]), float(plat[j])])
    frames = []
    for y in YEARS:
        ny = 0
        for m in range(1, 13):                       # one month per request (a full year exceeds the GEE user memory limit)
            t0, t1 = f"{y}-{m:02d}-01", (f"{y}-{m + 1:02d}-01" if m < 12 else f"{y + 1}-01-01")
            col = ee.ImageCollection(ASSET).filterDate(t0, t1).select(["sm_surface", "sm_rootzone"])
            rows = None
            for attempt in range(6):
                try:
                    rows = col.getRegion(pt, 9000).getInfo(); break
                except Exception as e:
                    print(f"  {k} {t0} attempt {attempt + 1} failed: {str(e)[:80]}", flush=True); time.sleep(10 * (attempt + 1))
            if rows is None:
                raise RuntimeError(f"{k} {t0}: all GEE attempts failed")
            df = pd.DataFrame(rows[1:], columns=rows[0])
            df["t"] = pd.to_datetime(df["time"], unit="ms")
            frames.append(df[["t", "sm_surface", "sm_rootzone"]]); ny += len(df)
        dy = pd.concat(frames[-12:])
        print(f"  {k} {y}: {ny} steps, surface median {dy.sm_surface.median():.3f}, root median {dy.sm_rootzone.median():.3f}", flush=True)
    d = pd.concat(frames).drop_duplicates("t").set_index("t").sort_index()
    if times is None:
        times = d.index
    d = d.reindex(times)
    SURF.append(d["sm_surface"].to_numpy(float)); ROOT.append(d["sm_rootzone"].to_numpy(float))

idx = np.array(list(sel.values())); bank = pbank[idx]
np.savez(f"{N}/noah/data/l4_points.npz", time=times.values.astype("datetime64[s]"), sm_surface=np.stack(SURF, 1),
         sm_rootzone=np.stack(ROOT, 1), name=np.array(list(sel.keys())), pix=idx, lat=plat[idx], lon=plon[idx],
         bank=bank, igbp=LC[bank], poros=B["poros"][bank], bee=B["bee"][bank], psis=B["psis"][bank], asset=ASSET)
print(f"saved noah/data/l4_points.npz: {len(times)} times x {len(idx)} pixels; NaN surface {np.isnan(np.stack(SURF)).mean():.2%}")
