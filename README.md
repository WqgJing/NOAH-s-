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

## Build environment (not built yet)
- Modules available: gcc, netcdf-fortran/4.6.1-mva2-hdf5-1.14, netcdf-c/4.9.2, mvapich2, openmpi.
- Build through sbatch, consistent with the other SLURM jobs.

## Planned modifications (none done yet)
1. BTRAN -> Laio beta(theta_L4_rootzone; s*, theta_w), bypassing the layer sum.
2. Prescribe soil moisture from L4 each step (surface for the top layer, root zone below); skip Noah-MP's own update.
3. Soil parameters from L4's Campbell curve (poros, bee, psis) per pixel.
4. OptDynamicVeg = 8 with observed LAI (Yuan v6.1 / MODIS).
5. ERA5 hourly forcing (tower forcing first at flux sites); TemperatureRadSfc output at the MODIS view time.
6. Nuisance parameters: LST offset, conductance scale, soil-resistance scale, roughness scale.
