"""MODIS daily LST at the Noah-MP test pixels (noah/data/l4_points.npz: BV, VR, KG), 2018-2024, for comparing with Noah-MP TRAD.
Products: MYD11A1 (Aqua, ~13:30 / 01:30 local solar) and MOD11A1 (Terra, ~10:30 / 22:30), collection 061, 1 km.
QC (same rule as sbi_train/vegsbi/amazon_lst_fetch.py 'err2'): LST produced (bits 0-1 <= 01) and LST error <= 2 K (bits 6-7 <= 01).
Per pixel and day: mean of valid 1 km LST over the 0.25 deg cell centred on the pixel (lat +- 0.125, lon +- 0.125), share of the cell's
1 km pixels that are valid, and mean view time (local solar hour) and view angle (deg) over the valid pixels; day and night separately.
Out: noah/data/lst_points.npz  (date [n_days], arrays [sat, n_days, n_pix] for keys ld, ln, fd, fn, vtd, vtn, vad, van; sat = aqua, terra)."""
import json, time
import numpy as np, pandas as pd, ee

N = "/storage/project/r-rbras6-0/jjason3/proposal_drought"
KEY = "/storage/project/r-rbras6-0/jjason3/drought/gldasnsmap-dee8cada466e.json"
SATS = {"aqua": "MODIS/061/MYD11A1", "terra": "MODIS/061/MOD11A1"}
YEARS = range(2018, 2025)
L = np.load(f"{N}/noah/data/l4_points.npz", allow_pickle=True); names, lat, lon = L["name"], L["lat"], L["lon"]
info = json.load(open(KEY))
ee.Initialize(ee.ServiceAccountCredentials(info["client_email"], KEY), project=info["project_id"],
              opt_url="https://earthengine-highvolume.googleapis.com")


def qc_ok(q):
    return q.bitwiseAnd(3).lte(1).And(q.rightShift(6).bitwiseAnd(3).lte(1))


BANDS = {}                                                     # actual view-angle band names, read from the collection


def layers(i):
    md, mn = qc_ok(i.select("QC_Day")), qc_ok(i.select("QC_Night"))
    ld = i.select("LST_Day_1km").multiply(0.02).updateMask(md)                 # K
    ln = i.select("LST_Night_1km").multiply(0.02).updateMask(mn)
    vtd = i.select("Day_view_time").multiply(0.1).updateMask(md); vtn = i.select("Night_view_time").multiply(0.1).updateMask(mn)
    vad = i.select(BANDS["vad"]).subtract(65).updateMask(md); van = i.select(BANDS["van"]).subtract(65).updateMask(mn)
    vd = md.unmask(0).rename("fd"); vn = mn.unmask(0).rename("fn")
    return ee.Image.cat([ld.rename("ld"), ln.rename("ln"), vtd.rename("vtd"), vtn.rename("vtn"), vad.rename("vad"), van.rename("van"), vd, vn])


KEYS = ["ld", "ln", "fd", "fn", "vtd", "vtn", "vad", "van"]
dates = pd.date_range(f"{YEARS[0]}-01-01", f"{YEARS[-1]}-12-31", freq="D"); nd = len(dates); npx = len(names)
OUT = {k: np.full((len(SATS), nd, npx), np.nan, np.float32) for k in KEYS}
for s_i, (sat, cid) in enumerate(SATS.items()):
    bn = ee.ImageCollection(cid).filterDate("2019-06-01", "2019-06-02").first().bandNames().getInfo()
    BANDS["vad"] = [b for b in bn if b.startswith("Day_view_ang")][0]
    BANDS["van"] = [b for b in bn if b.startswith("Night_view_ang")][0]
    print(f"{sat}: bands {bn}; view-angle bands {BANDS}", flush=True)
    for p in range(npx):
        cell = ee.Geometry.Rectangle([float(lon[p]) - 0.125, float(lat[p]) - 0.125, float(lon[p]) + 0.125, float(lat[p]) + 0.125])
        for y in YEARS:
            col = ee.ImageCollection(cid).filterDate(f"{y}-01-01", f"{y + 1}-01-01")
            fc = col.map(lambda i: ee.Feature(None, layers(i).reduceRegion(ee.Reducer.mean(), cell, 1000, maxPixels=1e8))
                         .set("date", i.date().format("YYYY-MM-dd")))
            res = None
            for attempt in range(6):
                try:
                    res = fc.getInfo(); break
                except Exception as e:
                    print(f"  {sat} {names[p]} {y} attempt {attempt + 1}: {str(e)[:80]}", flush=True); time.sleep(10 * (attempt + 1))
            if res is None:
                raise RuntimeError(f"{sat} {names[p]} {y}: all GEE attempts failed")
            for f in res["features"]:
                pr = f["properties"]; d = dates.get_loc(pd.Timestamp(pr["date"]))
                for k in KEYS:
                    v = pr.get(k)
                    if v is not None: OUT[k][s_i, d, p] = v
            m = (dates.year == y)
            print(f"{sat} {names[p]} {y}: days with valid day LST {np.isfinite(OUT['ld'][s_i, m, p]).sum()}, "
                  f"median clear share {np.nanmedian(OUT['fd'][s_i, m, p]):.2f}, median day LST {np.nanmedian(OUT['ld'][s_i, m, p]) - 273.15:.1f} C, "
                  f"median view time {np.nanmedian(OUT['vtd'][s_i, m, p]):.2f} h", flush=True)
np.savez(f"{N}/noah/data/lst_points.npz", date=dates.values.astype("datetime64[D]"), sat=np.array(list(SATS)), name=names, **OUT)
print("saved noah/data/lst_points.npz")
