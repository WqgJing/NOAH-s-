"""Compare the Python port (noah/pyport) with the Fortran reference runs runs/pyref_<NAME> (hpc/pyref_pixel.slurm).

Runs the Python model with the run directory's own inputs (namelist.hrldas, NoahmpTable.TBL, forcing/hrldas_setup.nc,
forcing/soil_spatial.nc, forcing/forcing_arrays.npz) and compares the half-hourly output with out/*.LDASOUT_DOMAIN1:
per variable max |diff|, RMS diff, the first time |diff| exceeds the per-step tolerance, and PASS/FAIL on the targets
(TRAD, TV, TG RMS <= 0.05 K; FSA, FIRA, LH, HFX, GRDFLX RMS <= 1 W m-2; BTRAN max <= 1e-4).
The first output record (ITIME = 0, before any physics) is skipped. ALBEDO is compared only where both are defined.
Writes the Python output to runs/pyref_<NAME>/pyport_out.npz.

Usage (SLURM only): sbatch --export=ALL,SCRIPT=pyport_compare.py hpc/plot.slurm
Environment: PYPORT_NAMES (default "VR KG"), PYPORT_NSTEPS (default: whole run), PYPORT_DETAIL (steps printed, 6).
"""
import glob
import os
import sys
import time

import numpy as np
import netCDF4

N = "/storage/project/r-rbras6-0/jjason3/proposal_drought/noah"
sys.path.insert(0, N)
from pyport.driver import PointModel  # noqa: E402

TEMPS = ["TRAD", "TV", "TG", "TGV", "TGB", "TAH", "T2MV", "T2MB", "SOIL_T1", "SOIL_T2", "SOIL_T3", "SOIL_T4"]
FLUXES = ["FSA", "FIRA", "SAV", "SAG", "LH", "HFX", "GRDFLX", "TR", "EVC", "IRC", "SHC", "IRG", "SHG", "EVG", "GHV",
          "IRB", "SHB", "EVB", "GHB", "APAR", "CANHS"]
OTHER_TOL = {"BTRAN": 1e-4, "COSZ": 1e-5, "SWFORC": 0.01, "LWFORC": 0.01, "RAINRATE": 1e-4, "ALBEDO": 1e-4,
             "EMISS": 1e-5, "FVEG": 1e-5, "LAI": 1e-5, "SAI": 1e-5, "EAH": 1.0, "Q2MV": 1e-5, "Q2MB": 1e-5,
             "CANLIQ": 1e-3, "FWET": 1e-3, "ECAN": 1e-6, "ETRAN": 1e-6, "EDIR": 1e-6, "CM": 1e-4, "CH": 1e-4,
             "CHV": 1e-4, "CHB": 1e-4, "CHLEAF": 1e-4, "CHUC": 1e-4, "CHV2": 1e-4, "CHB2": 1e-4, "PSN": 0.01,
             "SOIL_M1": 1e-4, "SOIL_M2": 1e-4, "SOIL_W1": 1e-4, "SOIL_W2": 1e-4}
TARGETS = {"TRAD": ("rms", 0.05), "TV": ("rms", 0.05), "TG": ("rms", 0.05), "FSA": ("rms", 1.0),
           "FIRA": ("rms", 1.0), "LH": ("rms", 1.0), "HFX": ("rms", 1.0), "GRDFLX": ("rms", 1.0),
           "BTRAN": ("max", 1e-4)}
REPORT = ["TRAD", "TV", "TG", "TGB", "TGV", "FSA", "FIRA", "LH", "HFX", "GRDFLX", "BTRAN", "SOIL_T1"]


def tolerance(name):
    if name in TEMPS:
        return 0.05
    if name in FLUXES:
        return 1.0
    return OTHER_TOL.get(name, 1e-4)


def read_ldasout(outdir, names):
    files = sorted(glob.glob(os.path.join(outdir, "*.LDASOUT_DOMAIN1")))
    times, data = [], {k: [] for k in names}
    for f in files:
        with netCDF4.Dataset(f) as nc:
            t = netCDF4.chartostring(nc["Times"][:])
            times += [np.datetime64(str(s).replace("_", "T"), "s") for s in np.atleast_1d(t)]
            for k in names:
                if k.startswith("SOIL_") and k[-1].isdigit():
                    a = nc[k[:-1]][:, 0, int(k[-1]) - 1, 0]
                elif k in nc.variables:
                    a = nc[k][:, 0, 0]
                else:
                    continue
                data[k].append(np.ma.filled(np.ma.asarray(a, dtype=np.float64), np.nan))
    data = {k: np.concatenate(v) for k, v in data.items() if v}
    return np.array(times), data


def compare(name, nsteps, ndetail):
    run = os.path.join(N, "runs", f"pyref_{name}")
    print(f"\n==================== {name}: {run}", flush=True)
    m = PointModel(run)
    p = m.p
    print(f"VegType {p.VegType} (MODIS), NROOT {p.NumSoilLayerRoot}, HVT {p.HeightCanopyTop}, Z0MVT {p.RoughLenMomVeg}, "
          f"VCMX25 {p.CarboxylRateMax25C}, s* {p.SoilMoistureLaioStar}, s_w {p.SoilMoistureLaioWilt}; "
          f"smcmax {p.SoilMoistureSat[0]:.4f} bexp {p.SoilExpCoeffB[0]:.3f} psisat {p.SoilMatPotentialSat[0]:.4f} "
          f"quartz {p.SoilQuartzFrac[0]:.3f}; start {m.start}, NTIME {m.ntime}", flush=True)
    nsteps = m.ntime if nsteps is None else min(nsteps, m.ntime)
    t0 = time.time()
    tp, py = m.run(nsteps)
    print(f"Python run: {nsteps} steps in {time.time() - t0:.1f} s", flush=True)
    b = py["BTRAN"]
    print(f"coverage: BTRAN < 1 in {np.mean(b < 1):.1%} of steps (min {b.min():.3g}, at the 1e-6 floor "
          f"{np.mean(b <= 1e-6):.1%}); LAI {py['LAI'].min():.3f}..{py['LAI'].max():.3f}; FVEG "
          f"{py['FVEG'].min():.3f}..{py['FVEG'].max():.3f}; rain steps {np.mean(py['RAINRATE'] > 0):.1%}; "
          f"wet canopy {np.mean(py['CANLIQ'] > 0):.1%}; night (COSZ <= 0) {np.mean(py['COSZ'] <= 0):.1%}; "
          f"min SOIL_T1 {py['SOIL_T1'].min():.2f} K", flush=True)
    np.savez(os.path.join(run, "pyport_out.npz"), time=tp, **py)

    tf, fo = read_ldasout(os.path.join(run, "out"), list(py.keys()))
    idx = {t: i for i, t in enumerate(tf)}
    missing = [t for t in tp if t not in idx]
    if missing:
        raise RuntimeError(f"{len(missing)} Python times not in LDASOUT, first {missing[0]}")
    sel = np.array([idx[t] for t in tp])

    print(f"\nfirst {ndetail} steps (Fortran / Python / diff):")
    for k in ["COSZ", "SWFORC", "FSA", "FIRA", "TRAD", "TV", "TG", "TAH", "EAH", "LH", "HFX", "GRDFLX", "BTRAN",
              "SOIL_T1"]:
        if k not in fo:
            continue
        f, q = fo[k][sel][:ndetail], py[k][:ndetail]
        print(f"  {k:8s} F " + " ".join(f"{x:11.5f}" for x in f))
        print(f"  {'':8s} P " + " ".join(f"{x:11.5f}" for x in q))
        print(f"  {'':8s} d " + " ".join(f"{x:11.2e}" for x in q - f))

    print(f"\n{'variable':10s} {'max|diff|':>11s} {'RMS diff':>11s} {'mean F':>11s} {'tol/step':>9s} "
          f"{'n>tol':>6s}  first time > tol")
    stats = {}
    for k in py:
        if k not in fo:
            continue
        f, q = fo[k][sel], py[k]
        ok = np.isfinite(f) & np.isfinite(q)
        if k == "ALBEDO":
            ok &= (f > -9998.0) & (q > -9998.0)
        d = q[ok] - f[ok]
        if d.size == 0:
            continue
        tol = tolerance(k)
        with np.errstate(invalid="ignore"):
            over = np.flatnonzero(ok & (np.abs(q - f) > tol))
        first = str(tp[over[0]]) if over.size else "-"
        mx, rms = float(np.max(np.abs(d))), float(np.sqrt(np.mean(d ** 2)))
        stats[k] = (mx, rms)
        print(f"{k:10s} {mx:11.4g} {rms:11.4g} {np.mean(f[ok]):11.4g} {tol:9.2g} {over.size:6d}  {first}")

    print(f"\n{name} report (max |diff|, RMS diff):")
    for k in REPORT:
        if k in stats:
            print(f"  {k:8s} max {stats[k][0]:10.4g}   RMS {stats[k][1]:10.4g}")
    allpass = True
    for k, (kind, lim) in TARGETS.items():
        val = stats[k][1] if kind == "rms" else stats[k][0]
        ok = val <= lim
        allpass &= ok
        print(f"  target {k:7s} {kind} {val:10.4g} <= {lim:g}: {'PASS' if ok else 'FAIL'}")
    print(f"{name}: {'PASS' if allpass else 'FAIL'}", flush=True)
    return allpass


if __name__ == "__main__":
    names = os.environ.get("PYPORT_NAMES", "VR KG").split()
    nsteps = int(os.environ["PYPORT_NSTEPS"]) if os.environ.get("PYPORT_NSTEPS") else None
    ndetail = int(os.environ.get("PYPORT_DETAIL", "6"))
    results = {n: compare(n, nsteps, ndetail) for n in names}
    print("\nSUMMARY: " + ", ".join(f"{n} {'PASS' if ok else 'FAIL'}" for n, ok in results.items()))
