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
- After each change: commit, then push to `git@github.com:WqgJing/NOAH-s-.git` (main).
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
- Cost: hrldas.exe takes 63 s per point-year. The driver's own timer reports SFLX (the land-model physics) as 0.00 s, so the time
  is file I/O (17,520 half-hourly forcing files plus outputs and restarts), not physics. For fitting s*, forcing must be kept in
  memory or the vector mode used.

## Build environment
- Modules available: gcc, netcdf-fortran/4.6.1-mva2-hdf5-1.14, netcdf-c/4.9.2, mvapich2, openmpi.
- Build through sbatch, consistent with the other SLURM jobs.

## Planned modifications (none done yet)
1. BTRAN -> Laio beta(theta_L4_rootzone; s*, theta_w), bypassing the layer sum.
2. Prescribe soil moisture from L4 each step (surface for the top layer, root zone below); skip Noah-MP's own update.
3. Soil parameters from L4's Campbell curve (poros, bee, psis) per pixel.
4. OptDynamicVeg = 8 with observed LAI (Yuan v6.1 / MODIS).
5. ERA5 hourly forcing (tower forcing first at flux sites); TemperatureRadSfc output at the MODIS view time.
6. Nuisance parameters: LST offset, conductance scale, soil-resistance scale, roughness scale.
