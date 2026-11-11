# pyport: Python port of our Noah-MP point configuration

This package ports the parts of Noah-MP that our point configuration uses. The source is HRLDAS v5.2.1 with Noah-MP and
`patches/noahmp_sstar.patch`, as built in `build/hrldas`. Everything runs in float64 NumPy for a single point. Each
function's docstring names the Fortran routine and file it ports, and variable names follow the Fortran (for example
`v.TemperatureCanopy`). The Fortran build is single precision (`kind_noahmp = 4`). To match its inputs exactly, the
port rounds all inputs, table parameters and physical constants to float32 and then computes in float64.

```python
import sys; sys.path.insert(0, "/storage/project/r-rbras6-0/jjason3/proposal_drought/noah")
from pyport import PointModel
m = PointModel("runs/pyref_VR")   # namelist.hrldas, NoahmpTable.TBL, forcing/{hrldas_setup.nc, soil_spatial.nc, forcing_arrays.npz}
times, out = m.run()              # half-hourly output with LDASOUT names (TRAD, TV, LH, ..., SOIL_T1..4)
```
A full run (122 days, 5856 steps) takes about 5 s. Python only runs through SLURM:
`sbatch --export=ALL,SCRIPT=pyport_compare.py hpc/plot.slurm`. Optional environment variables are
`PYPORT_NAMES="VR KG"`, `PYPORT_NSTEPS` and `PYPORT_DETAIL`.

## Configuration covered (checked against `runs/pyref_*/namelist.hrldas`; `driver.build_params` raises on any other)
These options are covered:
- DYNAMIC_VEG 8, CANOPY_STOMATAL_RESISTANCE 1 (Ball-Berry), BTR 4 (Laio), SURFACE_DRAG 1 (MOST)
- RADIATIVE_TRANSFER 3, SNOW_ALBEDO 1, SNOW_COVER 1, PCP_PARTITION 1
- SURFACE_RESISTANCE 4, TEMP_TIME_SCHEME 3, TBOT 2, SUPERCOOLED_WATER 1, SOIL_DATA 4
- PRESCRIBE_SOIL_MOISTURE 1; CROP / IRRIGATION / TILE / WETLAND 0
- NSOIL 4 (0.1/0.3/0.6/1.0 m), NOAH_TIMESTEP 1800 s, SOIL_TIMESTEP 0 (= model step), FORCING_TIMESTEP 3600 s
- ZLVL 10 m, MODIS (MODIFIED_IGBP_MODIS_NOAH) table; the `[noah/ s*]` nuisance scales are supported

## Ported routines
| Python (module.function) | Fortran |
|---|---|
| `driver.PointModel.init_state` | `hrldas/IO_code`: READINIT_HRLDAS, READVEG_HRLDAS (setup file). `drivers/hrldas`: NoahmpInitMain (cold start: TV=TG=TAH=TSK, EAH=2000, SMOIS cap, frozen-soil SH2O), NoahmpSnowInitMain (no snow) |
| `driver.HourlyForcing.at` | `module_hrldas_netcdf_io.F`: READFORC_HRLDAS, interpolate_inputstruct, copyfrom_inputstruct. T, Q, U, V, P, LW, SW, SMSURF, SMROOT are linear in time; RAINRATE and LAI are held from the earlier file |
| `driver.calc_declin` | `module_NoahMP_hrldas_driver.F`: CALC_DECLIN (cos SZA at the step's end time OLDDATE) |
| `driver.PointModel.step / read_forcing / prescribe_soil_moisture` | land_driver_exe: q to mixing ratio, RAINBL, first-step EAH/TAH/CH/CM guess, `[noah/ s*]` soil-moisture prescription (ice kept) |
| `driver.PointModel.column_vars / store_vars`, `build_params` | NoahmpDriverMain, ConfigVarInTransfer, ForcingVarInTransfer, EnergyVarInTransfer, WaterVarInTransfer, BiochemVarInTransfer, and the matching OutTransfer modules |
| `tables.*` | NoahmpReadTable (MODIS, rad, global, general, STAS soil blocks), NoahmpReadNamelist (defaults), READ_3D_SOIL |
| `main.noahmp_main` | NoahmpMain |
| `atmos.process_atmos_forcing`, `general_init`, `phenology_main` | AtmosForcingMod, GeneralInitMod, PhenologyMainMod |
| `canopy_water.canopy_water_intercept`, `precipitation_heat_advect`, `canopy_hydrology` | CanopyWaterInterceptMod, PrecipitationHeatAdvectMod, CanopyHydrologyMod |
| `main.energy_main` | EnergyMainMod (tile combination, TRAD, soil-temperature call, phase change, albedo) |
| `radiation.*` | SnowCoverGroundNiu07, SnowAgingBats, SnowAlbedoBats, GroundAlbedo, CanopyRadiationTwoStream, SurfaceAlbedo, SurfaceRadiation, SurfaceEmissivity |
| `ground.*` | GroundRoughnessProperty, SoilThermalProperty, GroundThermalProperty, SoilWaterTranspiration (option 4), ResistanceGroundEvaporation (options 1/4), PsychrometricVariable |
| `turbulence.*` | VaporPressureSaturation, ResistanceAboveCanopyMOST, ResistanceLeafToGround, ResistanceCanopyStomataBallBerry, SurfaceEnergyFluxVegetated, ResistanceBareGroundMOST, SurfaceEnergyFluxBareGround |
| `soil_thermal.*` | SoilSnowThermalDiffusion, MatrixSolverTriDiagonal, SoilSnowTemperatureSolver, SoilSnowTemperatureMain, SoilWaterSupercoolNiu06, SoilSnowWaterPhaseChange (soil layers) |
| `main.water_main_subset` | WaterMainMod: ground evaporation/dew split, CanopyHydrology, the frozen-ground top-layer ice update, and the no-snow checks of SnowWaterMain |

## Not ported, and why
- **Soil water solver**: SoilWaterMain (Richards equation, infiltration, runoff, drainage), groundwater and water-table
  routines. The prescription overwrites total soil water before every step, and these routines leave soil ice
  unchanged; ice is the only part of the soil water the prescription keeps. The end-of-step SOIL_M / SOIL_W from
  Python are therefore the values *before* the solver. They differ from LDASOUT (KG layer 1 by up to 0.025 m3/m3 after
  rain), but this has no effect on the next step.
- **Snow processes**: snow layers, compaction, snow hydrology, SNICAR and snow-on-ground temperature resets. The
  window is snow-free. The code raises `NotImplementedError` if snowfall reaches the ground, if there is initial
  snow, or if snow depth exceeds 0.05 m. BATS snow aging and albedo and the Niu07 snow-cover fraction are ported
  because they run every step, and they evaluate to zero snow.
- **Other land types and options**: glacier, lake, urban, crop model and irrigation. (2026-11-10: cropland now runs when IRRIGATION_OPTION = 0, verified at BV May-Sep 2019, TRAD RMS 3.8e-4 K; previously the code raised for cropland, because
  IrrigationPrepare can trigger irrigation there even with IRRIGATION_OPTION 0. Also excluded: wetland, tile
  drainage, dynamic-vegetation biochemistry, and every alternative option (Jarvis, Chen97, other BTR / RSURF /
  radiation / albedo schemes).
- **Balance checks and diagnostics**: BalanceWaterInit/Check and BalanceEnergyCheck are diagnostic only.
  IrrigationPrepare is a no-op for non-cropland.

## Verification (`scripts/pyport_compare.py` against Fortran LDASOUT, 2018-06-01 to 2018-10-01, 5856 half-hourly steps)
The Python model reads the reference run's own inputs: the same setup file (January-based initial TSLB/SMOIS/TSK),
the same forcing and the same namelist. In the first steps it agrees with the Fortran to single-precision roundoff:
TRAD differs by about 1e-5 K and LH by about 1e-4 W m-2. Over the whole window (max |diff| / RMS diff):

| variable | VR | KG | VRs (s* 0.25, s_w 0.10) | KGsn (s* 0.25, s_w 0.10, scales 1.5/0.7/1.3) |
|---|---|---|---|---|
| TRAD [K] | 0.0046 / 1.4e-4 | 0.0071 / 1.4e-4 | 0.0036 / 9.7e-5 | 0.0048 / 1.1e-4 |
| TV [K] | 0.0080 / 2.2e-4 | 0.0079 / 1.8e-4 | 0.0095 / 2.7e-4 | 0.0068 / 1.6e-4 |
| TG [K] | 0.0044 / 1.4e-4 | 0.0059 / 1.2e-4 | 0.0037 / 1.1e-4 | 0.0043 / 9.5e-5 |
| TGB [K] | 0.0041 / 6.3e-5 | 0.00039 / 3.1e-5 | 0.0042 / 6.3e-5 | 0.00023 / 2.2e-5 |
| TGV [K] | 0.0093 / 2.4e-4 | 0.017 / 3.6e-4 | 0.0060 / 1.8e-4 | 0.016 / 3.3e-4 |
| FSA [W m-2] | 0.0016 / 3.5e-5 | 0.00022 / 2.3e-5 | 0.0016 / 3.5e-5 | 0.00022 / 2.3e-5 |
| FIRA [W m-2] | 0.030 / 8.3e-4 | 0.052 / 9.6e-4 | 0.020 / 5.6e-4 | 0.036 / 7.6e-4 |
| LH [W m-2] | 0.022 / 3.7e-4 | 0.010 / 2.1e-4 | 0.038 / 6.7e-4 | 0.027 / 5.2e-4 |
| HFX [W m-2] | 0.084 / 0.0022 | 0.16 / 0.0030 | 0.046 / 0.0016 | 0.089 / 0.0019 |
| GRDFLX [W m-2] | 0.055 / 0.0017 | 0.11 / 0.0020 | 0.043 / 0.0014 | 0.053 / 0.0012 |
| BTRAN [-] | 0 / 0 | 0 / 0 | 5.7e-8 / 2.5e-8 | 3.4e-8 / 1.4e-8 |
| SOIL_T layer 1 [K] | 6.3e-4 / 5.7e-5 | 9.0e-4 / 5.7e-5 | 5.3e-4 / 5.6e-5 | 6.0e-4 / 4.1e-5 |

All targets pass with large margins at all four runs: TRAD, TV and TG RMS below 0.05 K; FSA, FIRA, LH, HFX and
GRDFLX RMS below 1 W m-2; BTRAN max below 1e-4. Across the 21 component energy terms (SAV, SAG, IRC, SHC, EVC, TR,
IRG, SHG, EVG, GHV, IRB, SHB, EVB, GHB, ...) the largest single-step difference is 0.41 W m-2 (SHG, KG). The rare
per-step exceedances of the diagnostic tolerances are float32 effects at thresholds: T2MV / Q2MV when the 2-m
exchange coefficient is near 1e-5, FWET at the 1e-6 mm canopy-water cutoff, and the CH* conductances.

**What these runs cover.** In the stock runs BTRAN = 1 at every step, because root-zone moisture stays above s*. The
variant runs VRs and KGsn therefore use the same inputs with s* = 0.25 and s_w = 0.10 (BTRAN 0.27–0.99 at VR and
0.085–0.99 at KG). KGsn also sets the nuisance scales to 1.5 / 0.7 / 1.3. These variants are made with
`hpc/pyref_variant.slurm` and `scripts/make_pyref_variant.py`. Across the runs LAI is 0.35–1.51, rain falls in 20–37 %
of steps and the canopy is wet in 2–12 %.

**Not exercised by any reference run** (ported, not verified):
- frozen soil and the phase change with ice (minimum soil temperature in the window is 283 K)
- the frozen-soil SH2O initialization
- the bare-only branch (LAI < 0.05)
- the BTRAN floor of 1e-6
- canopy snow / ice
- TG at or below freezing

## Sensitivity hook (Python only, not in the Fortran)
`m.p.BtranOverride = b` fixes BTRAN at b (ignores soil moisture); default None. s* can be changed with
`m.p.SoilMoistureLaioStar = s` after building the model. With the default (None) the output is bit-identical to the verified
run (checked in scripts/sens_trad_beta.py).
