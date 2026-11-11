# noah/ — Noah-MP workspace for the LST-based s* estimate

Purpose: forward model of LST with Laio's water-stress function on the canopy, to infer the root-zone threshold s*
from MODIS LST. The design discussion is in `../tasks/todo.md` (2026-11-10 sections).

## Source (fetched 2026-11-10, untouched upstream code)
- `src/hrldas/`: NCAR HRLDAS (offline driver), tag **v5.2.1**, from https://github.com/NCAR/hrldas
- `src/hrldas/noahmp/`: NCAR Noah-MP submodule, commit **17751dcd7a2442a3c57138f530f523daf08dc1c5**, from https://github.com/NCAR/noahmp
- Our own scripts and modifications will go in a separate folder (e.g. `noah/scripts/`, `noah/patches/`), not inside `src/`.
- `src/hrldas` is a git submodule pinned to v5.2.1 (commit 8e07c72), so this repo stores only the pointer, not NCAR's code.
  Clone with `git clone --recurse-submodules <this repo>`; an existing clone needs `git submodule update --init --recursive`.

## Code facts checked against the source (noahmp/src)
- **BTRAN** = `SoilTranspFacAcc` (SoilWaterTranspirationMod.F90).
  - Per root layer `SoilWetFac`:
    - opt 1 Noah: (theta - theta_wilt)/(theta_fc - theta_wilt), clipped to [0,1], i.e. Laio's linear shape with table breakpoints;
    - opt 2 CLM and opt 3 SSiB: functions of matric potential.
  - The accumulated value is a THICKNESS-weighted sum over the root-zone layers, not root-fraction weighted.
- **Where BTRAN acts:**
  - Ball-Berry: multiplies Vcmax (ResistanceCanopyStomataBallBerryMod.F90:118).
  - Jarvis: rs = rs_min / (f_solar f_temp f_vpd BTRAN) (ResistanceCanopyStomataJarvisMod.F90:93).
  - Dynamic vegetation: SoilWaterStress = 1 - BTRAN (BiochemNatureVegMainMod.F90:95), so dynamic veg must be OFF.
- **OptDynamicVeg = 8**: input LAI, VegFrac = 1 - exp(-0.52 (LAI + SAI)), with SAI = max(0.05, 0.1 LAI) (PhenologyMainMod.F90:94-138).
  VegFrac has a floor of 0.05.
- **Soil evaporation resistance** (OptGroundResistanceEvap): 1 Sakaguchi & Zeng 2009, 2 Sellers 1992, 3 adjusted Sellers, 4 opt 1 + snow.
- **Radiative temperature output**: `TemperatureRadSfc` (EnergyVarType.F90:207).
- **Campbell soil parameters**: SoilMoistureSat, SoilExpCoeffB, SoilMatPotentialSat (match L4 poros, bee, psis).
- **Examples**: hrldas/HRLDAS_forcing/run/examples/{single_point (Bondville), ERA5 (ERA5-Land forcing prep)}.

## Workflow
- After each change: commit. Push to `git@github.com:WqgJing/NOAH-s-.git` (main) for real changes (model code, setup,
  results). Small checks and plot tweaks are committed locally only, not pushed.
- Upstream code under `src/` is never edited in place. Our changes go in `scripts/`, `hpc/` and (later) `patches/`.
- Figures: line plots only, minimal and clean. No titles or explanatory text; axis labels with units, and a legend only if
  there is more than one line.
- Jobs run through `sbatch` (`hpc/*.slurm`, logs in `logs/`). Batch scripts must `source /etc/profile.d/z00_lmod_pace.sh`
  before `module load`.

## Reproduction of NCAR's single-point example (2026-11-10): PASSED
- Build: `sbatch hpc/build_hrldas.slurm` (serial gfortran 12.3, netCDF-Fortran 4.6.1; only the netCDF paths changed), 49 s.
- Run: `sbatch --export=ALL,CONFIG=shipped|tutorial hpc/run_bondville.slurm`. This is Bondville IL, 1998, 30-min forcing.
  Starts 1998-01-02 (the forcing starts then, as in the tutorial), KDAY 364, daily output files.
- Compare: `sbatch hpc/compare_bondville.slurm` writes `runs/figs/bondville_reproduce.png`. References were extracted from
  `src/hrldas/tutorial/Note1_Single_Point.ipynb`.
- Results:
  - LH on 20 Jun 1998 vs the tutorial's ncview screenshot (digitised): r 1.000, RMSE 1.9 W m-2, peak 428 vs 432 at 18:30 UTC.
    The shipped and tutorial options give the same result here.
  - Top-soil temperature on 2 Jan 00:00-07:30 vs the tutorial run log: max difference 0.10 K (tutorial options) and 0.27 K
    (shipped options). The remainder is most likely the tutorial's older release: its namelist lacks several v5.2 options.
  - Year: ET about 720 mm, Jun-Aug mean LH 105 W m-2, TRAD 247-310 K. The first output step carries -9999 fill values (masked).
- Cost: hrldas.exe takes 63-72 s per point-year, of which the land-model physics is only 1.3-1.4 s (cumulative SFLX timer
  added in the patch; debug flags -g -fbounds-check). The rest is file I/O (17,520 half-hourly forcing files plus outputs).
  CORRECTION: an earlier note read the stock "SFLX: 0.00" as the year's physics time; that timer is reset every step.
  For fitting s*, forcing must be kept in memory or the vector mode used.

## Build environment
- Modules available: gcc, netcdf-fortran/4.6.1-mva2-hdf5-1.14, netcdf-c/4.9.2, mvapich2, openmpi.
- Build through sbatch, consistent with the other SLURM jobs.

## Our modifications (patches/noahmp_sstar.patch, 15 files; every change is marked "[noah/ s*]" in the code)
Work flow: `scripts/prepare_build.sh` (build/hrldas = src/hrldas + patch) -> `sbatch hpc/build_patched.slurm`.
After editing build/hrldas: `scripts/make_patch.sh` re-exports the patch. src/ is never edited.

| Namelist switch (&NOAHLSM_OFFLINE) | Default (= stock) | When on |
|---|---|---|
| `BTR_OPTION = 4` | 1 | BTRAN = Laio beta = clip((theta_root - LAIO_SWILT)/(LAIO_SSTAR - LAIO_SWILT), 0, 1). theta_root = prescribed root zone (below), else the thickness-weighted root-layer mean liquid water. Uptake spread by layer thickness. |
| `LAIO_SSTAR`, `LAIO_SWILT` | -9999 (must be set for option 4) | s*, s_w [m3/m3]; checked 0 <= s_w < s* <= 1 |
| `PRESCRIBE_SOIL_MOISTURE = 1` | 0 | each step before the physics: layer 1 <- forcing SMSURF, layers 2..NSOIL <- SMROOT (existing ice kept) |
| `NUIS_GS_SCALE` | 1.0 | Vcmax25 x s (Ball-Berry), minimum stomatal resistance / s (Jarvis) |
| `NUIS_RSURF_SCALE` | 1.0 | soil evaporation resistance x s (soil formula, before snow mixing) |
| `NUIS_Z0_SCALE` | 1.0 | canopy roughness length, canopy top and bottom height x s |
| (always) | - | new output `BTRAN`; cumulative physics timer `SFLXtot` in the log |

Upstream bug fixed in the same patch (ConfigVarInTransferMod.F90): with `SOIL_DATA_OPTION = 4` the per-layer soil type stayed
undefined (-9999) while EnergyVarInTransfer / WaterVarInTransfer index the soil tables with it before the option-4 values
overwrite them, so the run aborts under -fbounds-check ("Index '-9999' ... quartz_table"). Now initialized from ISLTYP as option 3
does. Options 1-3 unchanged (regression test still bit-identical).

No code is needed for observed LAI (`DYNAMIC_VEG_OPTION = 8`: LAI from each forcing file, vegetation fraction from LAI) or for
L4 soil parameters (`SOIL_DATA_OPTION = 4` + `SPATIAL_FILENAME`).

Tests (`sbatch hpc/test_patched_bondville.slurm`, Bondville 1998; all PASS, 2026-11-10):
1. All switches off vs the stock run: bit-identical (36,135 variable-files, 0 differences).
2. `BTR_OPTION 4` without prescription: BTRAN = Laio(start-of-step root-zone liquid water) to 2.8e-7 (float32). Cropland has 3
   root layers in the MODIS table (0-1 m, the same depth as the L4 root zone).
3. Prescription: BTRAN = Laio(prescribed SMROOT) exactly (lag 0); SOIL_M follows SMSURF / SMROOT (median 0.0003 / 0.0001).
4. Nuisance directions (Jun-Aug daytime): conductance x2 -> transpiration +21 %, TRAD -0.39 K; soil resistance x2 -> soil
   evaporation -6 %; roughness x2 -> TRAD -0.06 K.

## Inputs for real pixels (scripts/make_pixel_forcing.py, hpc/run_pixel.slurm)
- ERA5 hourly (clima/forcing_global/era5_hourly; longwave msdwlwrf and precipitation mtpr added 2026-11-10 with
  clima/forcing_data/scripts/fetch_era5_noah.slurm), SW/LW centred on the hour, precipitation applied over the next hour.
- SMAP L4 SPL4SMGP/008 3-hourly surface + root zone at the test pixels (scripts/fetch_l4_points.py), Yuan v6.1 LAI.
- Soil: L4 Campbell (porosity, b, psi_s) in soil_spatial.nc; wilting -153 m and field capacity -3.36 m on that curve; quartz /
  Ks placeholders from the loam table class (theta is prescribed; quartz affects soil heat conduction, to revisit).
- Laio: s* = bank theta_star_ref, s_w = bank theta_w.
- MODIS Aqua/Terra daily LST at the same pixels (scripts/fetch_lst_points.py) for the TRAD comparison (scripts/compare_trad_lst.py).
- SMAP L4 root zone comes from the local daily files (sbi_train/runs/l4/l4_rootzone_<year>.npy, the map's product). They are
  identical to GEE 008 for 2018-2022 (scripts/check_l4_local.py). The surface layer is from GEE 008 (no local L4 surface).

## First forward runs (2018 spin-up, 2019 shown; hpc/run_pixel.slurm, scripts/plot_pixel_runs.py -> runs/figs/pixel_runs_2019.png)
Reference s* (bank theta_star_ref), not fitted:
- BV (cropland, s* 0.298): L4 root zone below s* on 25 % of days (Aug-Oct); daytime BTRAN down to 0.57.
- VR (s* 0.129) and KG (s* 0.077), grassland: L4 root zone never below s* (minimum 0.141, 0.114), so BTRAN = 1 all year and there
  is no s* signal in TRAD at this s*.
- TRAD range 2019: BV -35 to 41 C (Jan 2019 cold outbreak), VR -4 to 53 C, KG -7 to 60 C.
- Run time 5.5-7.5 min per pixel for 2 years; physics 2.5-3.9 s, the rest file I/O.

## Python port (pyport/, see pyport/README.md)
NumPy port of the routines this configuration uses (driver, radiation, two-stream, vegetated and bare energy balance, MOST,
Ball-Berry, Laio BTRAN, soil heat, phase change, canopy water). Snow, glacier, crop, irrigation, runoff/groundwater and the
soil-water solver (overwritten by the prescription) are excluded; the code raises if they would be needed.
Verified against Fortran reference runs (hpc/pyref_pixel.slurm, Jun-Sep 2018, VR and KG, plus variants with s* = 0.25 so that the
Laio branch is active; scripts/pyport_compare.py): TRAD RMS diff about 1e-4 K, fluxes about 1e-3 W m-2 RMS, BTRAN 3e-8, i.e.
single-precision round-off. About 5 s per 122-day run. Time series: runs/figs/pyport_vs_fortran.png (scripts/plot_pyport_compare.py).
Gaps (details in pyport/README.md):
- Ported but not exercised by any reference run: frozen soil / phase change with ice, the frozen-soil start-up formula, the
  no-vegetation (bare-only) branch, the BTRAN floor 1e-6, canopy snow / ice.
- Not ported, the code raises: snow on the ground (layers, compaction, snow hydrology); cropland (the irrigation trigger can fire
  even with irrigation off).
- Left out by design: soil-water solver, runoff, groundwater (overwritten by the L4 prescription); glacier, lake, urban, wetland,
  tile drainage, dynamic vegetation, alternative schemes.
- The full Fortran model has none of these limits. Next check: the port on the 2018-2019 VR / KG runs (winter frost nights).

## TRAD sensitivity to beta (scripts/sens_trad_beta.py -> runs/figs/sens_trad_beta.png; Python port, Jun-Sep 2018)
TRAD at the Aqua overpass (13:30 local), seasonal mean change vs beta = 1:

| | VR (FVEG ~0.5) | KG (FVEG ~0.2) |
|---|---|---|
| Transpiration at beta = 1 | 1.0 mm/d | 0.25 mm/d |
| Fixed beta 1 -> 0.05: TRAD | +0.45 K | +0.13 K |
| Fixed beta 1 -> 0.05: canopy TV | +0.92 K | +0.56 K |
| s* scan up to "always below s*": TRAD, mean (max) | +0.17 (0.46) K | +0.02 (0.07) K |

So in the dry season at sparse grassland the closure signal in TRAD is below MODIS LST error (1-2 K): little transpiration is
left to switch off, and only the canopy part (weighted by FVEG) warms. Signal should be larger where and when transpiration is
high (spring drydown, cropland, denser vegetation); not tested yet.

### Growing season 2019 (scripts/sens_beta_year_py.py -> runs/figs/sens_beta_year_py.png; Python port, 1 Mar - 31 Oct, cold start)
- The full-year attempt stopped in Jan 2018 at both pixels: the port refuses any snowfall (ERA5 had ~1e-9 mm/s). Until then it
  matched Fortran (TRAD RMS 8e-5 / 4e-5 K).
- dTRAD at 13:30 for fixed beta 0.05 vs 1 is proportional to the transpiration heat available at 13:30: about 1.4 K per
  100 W m-2 at VR (r near 1 in the scatter).
- VR: April mean +1.4 K (transpiration heat 100 W m-2, LAI 2.3, FVEG 0.72); daily peaks +2.6 K in late April (180 W m-2);
  Jul-Oct +0.3 to +0.5 K.
- KG: at most +0.2 K (monsoon, transpiration heat at most about 20 W m-2).
- So TRAD can only carry s* where midday transpiration is well above 100 W m-2 (green, dense vegetation): about 1.5 K per
  100 W m-2 here.

### Cropland, Bondville (BV) May-Sep 2019 (Python port; reference hpc/pyref_pixel.slurm START=2019-05-01 KDAY=153)
- Cropland enabled in the port when IRRIGATION_OPTION = 0: the irrigation trigger sets FlagIrri = .false. and the irrigation
  evaporation heat is reset to 0 every step, so the cropland branches do nothing. Port vs Fortran at BV: all targets PASS
  (TRAD RMS 3.8e-4 K, LH 0.008 W m-2, BTRAN 3e-8), with BTRAN < 1 on 38 % of steps.
- Fixed-beta sweep (runs/figs/sens_beta_year_py_BV.png), monthly mean dTRAD at 13:30:

  | Month | LAI | Transpiration heat 13:30 | beta 0.5 | beta 0.2 | beta 0.05 |
  |---|---|---|---|---|---|
  | Jul | 2.9 | 233 W m-2 | +1.0 K | +2.3 K | +3.0 K |
  | Aug | 3.6 | 278 W m-2 | +1.3 K | +2.7 K | +3.5 K |

  Daily peaks about +6 K (beta 0.05) and +3.7 K (beta 0.5) at 400 W m-2. The same slope as VR, about 1.4-1.5 K per
  100 W m-2. A green crop in midsummer is where TRAD is clearly sensitive to beta.
