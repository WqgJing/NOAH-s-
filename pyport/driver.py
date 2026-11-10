"""HRLDAS driver behaviour for one point: parameters, initialization, forcing time interpolation, solar zenith angle,
soil-moisture prescription, the per-step transfer into / out of the column model, and the time loop.

Ports (for our configuration):
  * hrldas/IO_code/module_NoahMP_hrldas_driver.F: land_driver_ini (cold start), land_driver_exe (time loop body:
    forcing unit conversions, EAH/TAH/CH/CM first-step guess, [noah/ s*] soil-moisture prescription), CALC_DECLIN
  * hrldas/IO_code/module_hrldas_netcdf_io.F: READFORC_HRLDAS / interpolate_inputstruct / copyfrom_inputstruct
  * noahmp/drivers/hrldas/NoahmpInitMainMod.F90 (NoahmpInitMain), NoahmpSnowInitMod.F90 (no snow)
  * noahmp/drivers/hrldas/NoahmpDriverMainMod.F90 and the *VarInTransfer / *VarOutTransfer modules
"""
import math
import os
from datetime import datetime, timedelta
from types import SimpleNamespace

import numpy as np

from .constants import f32
from . import tables
from .main import noahmp_main

F32 = np.float32


# ---------------------------------------------------------------------------------------------------------------------
# parameters and options
# ---------------------------------------------------------------------------------------------------------------------
def build_params(nl, tbl, setup, soil):
    """Options, domain constants and parameters as set by ConfigVarInTransfer, EnergyVarInTransfer,
    WaterVarInTransfer and BiochemVarInTransfer (constant in time for one point)."""
    p = SimpleNamespace()
    # options (NoahmpReadNamelistMod.F90 -> ConfigVarInTransferMod.F90)
    p.OptDynamicVeg = int(nl["DYNAMIC_VEG_OPTION"])
    p.OptStomataResistance = int(nl["CANOPY_STOMATAL_RESISTANCE_OPTION"])
    p.OptSoilWaterTranspiration = int(nl["BTR_OPTION"])
    p.OptSurfaceDrag = int(nl["SURFACE_DRAG_OPTION"])
    p.OptSoilSupercoolWater = int(nl["SUPERCOOLED_WATER_OPTION"])
    p.OptSoilPermeabilityFrozen = int(nl["FROZEN_SOIL_OPTION"])   # soil water only (not ported, no effect here)
    p.OptCanopyRadiationTransfer = int(nl["RADIATIVE_TRANSFER_OPTION"])
    p.OptSnowAlbedo = int(nl["SNOW_ALBEDO_OPTION"])
    p.OptSnowThermConduct = int(nl["SNOW_THERMAL_CONDUCTIVITY"])
    p.OptSnowCoverGround = int(nl["SNOW_COVER_OPTION"])
    p.OptRainSnowPartition = int(nl["PCP_PARTITION_OPTION"])
    p.OptSoilTemperatureBottom = int(nl["TBOT_OPTION"])
    p.OptSnowSoilTempTime = int(nl["TEMP_TIME_SCHEME_OPTION"])
    p.OptGroundResistanceEvap = int(nl["SURFACE_RESISTANCE_OPTION"])
    p.OptSoilProperty = int(nl["SOIL_DATA_OPTION"])
    p.OptCropModel = int(nl["CROP_OPTION"])
    p.OptIrrigation = int(nl["IRRIGATION_OPTION"])
    p.OptTileDrainage = int(nl["TILE_DRAINAGE_OPTION"])
    p.OptWetlandModel = int(nl["WETLAND_OPTION"])
    p.PrescribeSoilMoisture = int(nl["PRESCRIBE_SOIL_MOISTURE"])
    req = dict(OptDynamicVeg=8, OptStomataResistance=1, OptSoilWaterTranspiration=4, OptSurfaceDrag=1,
               OptCanopyRadiationTransfer=3, OptSnowAlbedo=1, OptGroundResistanceEvap=4, OptSnowSoilTempTime=3,
               OptSoilTemperatureBottom=2, OptSoilSupercoolWater=1, OptSoilProperty=4, PrescribeSoilMoisture=1,
               OptCropModel=0, OptIrrigation=0, OptTileDrainage=0, OptWetlandModel=0, OptRainSnowPartition=1,
               OptSnowCoverGround=1)
    bad = {k: (getattr(p, k), want) for k, want in req.items() if getattr(p, k) != want}
    if bad:
        raise NotImplementedError(f"options differ from the ported configuration (got, ported): {bad}")
    if int(nl.get("SF_URBAN_PHYSICS", 0)) != 0:
        raise NotImplementedError("urban physics is not ported")

    # domain
    p.NumSoilLayer = int(nl["NSOIL"])
    p.NumSwRadBand = 2
    p.MainTimeStep = f32(nl["NOAH_TIMESTEP"])
    soil_steps = max(int(round(nl["SOIL_TIMESTEP"] / p.MainTimeStep)), 1)
    if soil_steps != 1:
        raise NotImplementedError("SOIL_TIMESTEP different from NOAH_TIMESTEP is not ported")
    p.NumSoilTimeStep = soil_steps
    p.SoilTimeStep = p.MainTimeStep * soil_steps
    p.FlagSoilProcess = True
    dzs = nl["SOIL_THICK_INPUT"]
    zsoil = np.zeros(p.NumSoilLayer, dtype=F32)                    # ZSOIL computed in float32 as in the driver
    zsoil[0] = -F32(dzs[0])
    for k in range(1, p.NumSoilLayer):
        zsoil[k] = -F32(dzs[k]) + zsoil[k - 1]
    p.DepthSoilLayer = [float(x) for x in zsoil]
    p.RefHeightAboveSfc = float(F32(2.0) * F32(nl["ZLVL"]) * F32(0.5))      # DZ8W = 2 ZLVL; RefHeight = DZ8W * 0.5
    p.DepthSoilTempBottom = tbl["ZBOT"]
    p.VegType = int(setup["IVGTYP"])
    p.SoilTypeTable = int(setup["ISLTYP"])
    p.IndexWaterPoint, p.IndexBarrenPoint, p.IndexIcePoint = tbl["ISWATER"], tbl["ISBARREN"], tbl["ISICE"]
    p.FlagUrban = (p.VegType == tbl["ISURBAN"]) or (p.VegType > tbl["URBTYPE_BEG"])
    if p.FlagUrban or p.VegType in (p.IndexWaterPoint, p.IndexIcePoint):
        raise NotImplementedError("urban / water / glacier points are not ported")
    p.Latitude = setup["XLAT"]

    # energy parameters (EnergyVarInTransferMod.F90), nuisance scales applied as in the patch
    p.SoilHeatCapacity = tbl["CSOIL"]
    p.SnowAgeFacBats, p.SnowGrowVapFacBats = tbl["TAU0"], tbl["GRAIN_GROWTH"]
    p.SnowSootFacBats, p.SnowGrowFrzFacBats = tbl["DIRT_SOOT"], tbl["EXTRA_GROWTH"]
    p.SolarZenithAdjBats = tbl["BATS_COSZ"]
    p.FreshSnoAlbVisBats, p.FreshSnoAlbNirBats = tbl["BATS_VIS_NEW"], tbl["BATS_NIR_NEW"]
    p.SnoAgeFacDifVisBats, p.SnoAgeFacDifNirBats = tbl["BATS_VIS_AGE"], tbl["BATS_NIR_AGE"]
    p.SzaFacDirVisBats, p.SzaFacDirNirBats = tbl["BATS_VIS_DIR"], tbl["BATS_NIR_DIR"]
    p.UpscatterCoeffSnowDir, p.UpscatterCoeffSnowDif = tbl["BETADS"], tbl["BETAIS"]
    p.ZilitinkevichCoeff = tbl["CZIL"]
    p.EmissivitySnow = tbl["SNOW_EMIS"]
    p.EmissivitySoilLake = tbl["EG"]
    p.RoughLenMomSnow, p.RoughLenMomSoil = tbl["Z0SNO"], tbl["Z0SOIL"]
    p.ResistanceSoilExp, p.ResistanceSnowSfc = tbl["RSURF_EXP"], tbl["RSURF_SNOW"]
    p.ResistanceGrdEvapScale = nl["NUIS_RSURF_SCALE"]
    z0s = nl["NUIS_Z0_SCALE"]
    p.HeightCanopyTop = f32(tbl["HVT"] * z0s)
    p.HeightCanopyBot = f32(tbl["HVB"] * z0s)
    p.RoughLenMomVeg = f32(tbl["Z0MVT"] * z0s)
    p.TreeCrownRadius = tbl["RC"]
    p.CanopyWindExtFac = tbl["CWPVT"]
    p.CanopyOrientIndex = tbl["XL"]
    p.ConductanceLeafMin = tbl["BP"]
    p.Co2MmConst25C, p.O2MmConst25C = tbl["KC25"], tbl["KO25"]
    p.Co2MmConstQ10, p.O2MmConstQ10 = tbl["AKC"], tbl["AKO"]
    p.LeafDimLength = tbl["DLEAF"]
    p.HeatCapacCanFac = tbl["CBIOM"]
    p.ReflectanceLeaf, p.ReflectanceStem = tbl["RHOL"], tbl["RHOS"]
    p.TransmittanceLeaf, p.TransmittanceStem = tbl["TAUL"], tbl["TAUS"]
    p.AlbedoSoilSat, p.AlbedoSoilDry = tbl["ALBSAT"], tbl["ALBDRY"]      # soil colour 4 (ConfigVarInTransferMod.F90)
    p.ScatterCoeffSnow = tbl["OMEGAS"]
    p.CO2, p.O2 = tbl["CO2"], tbl["O2"]

    # water parameters (WaterVarInTransferMod.F90), SOIL_DATA_OPTION 4 = spatial soil file
    p.CanopyLiqHoldCap = tbl["CH2OP"]
    p.SnowfallDensityMax = tbl["SNOWDEN_MAX"]
    p.SnowMassFullCoverOld = tbl["SWEMX"]
    p.SnowMeltFac, p.SnowCoverFac = tbl["MFSNO"], tbl["SCFFAC"]
    p.NumSoilLayerRoot = tbl["NROOT"]
    p.SoilMoistureLaioStar = nl["LAIO_SSTAR"]
    p.SoilMoistureLaioWilt = nl["LAIO_SWILT"]
    if not (0.0 <= p.SoilMoistureLaioWilt < p.SoilMoistureLaioStar <= 1.0):
        raise ValueError("LAIO_SWILT / LAIO_SSTAR invalid")
    n = p.NumSoilLayer
    p.SoilExpCoeffB = list(soil["bexp"][:n])
    p.SoilMoistureWilt = list(soil["smcwlt"][:n])
    p.SoilMoistureFieldCap = list(soil["smcref"][:n])
    p.SoilMoistureSat = list(soil["smcmax"][:n])
    p.SoilMatPotentialSat = list(soil["psisat"][:n])
    p.SoilQuartzFrac = list(soil["quartz"][:n])

    # biochemistry parameters (BiochemVarInTransferMod.F90)
    p.NitrogenConcFoliageMax = tbl["FOLNMX"]
    p.QuantumEfficiency25C = tbl["QE25"]
    p.CarboxylRateMax25C = f32(tbl["VCMX25"] * nl["NUIS_GS_SCALE"])
    p.CarboxylRateMaxQ10 = tbl["AVCMX"]
    p.PhotosynPathC3 = tbl["C3PSN"]
    p.SlopeConductToPhotosyn = tbl["MP"]
    p.TemperatureMinPhotosyn = tbl["TMIN"]
    return p


# ---------------------------------------------------------------------------------------------------------------------
# forcing
# ---------------------------------------------------------------------------------------------------------------------
class HourlyForcing:
    """Forcing at the model time, as READFORC_HRLDAS returns it (module_hrldas_netcdf_io.F).

    At a forcing time: the values of that file. Between two forcing times A < t < B: T, Q, U, V, PSFC, LWDOWN,
    SWDOWN and [noah/ s*] SMSURF, SMROOT linear in time (weight of A = (B - t)/(B - A)); RAINRATE and LAI held at the
    values of A. All values are float32 as in the LDASIN files and the Fortran arrays."""
    LINEAR = ("T2D", "Q2D", "U2D", "V2D", "PSFC", "LWDOWN", "SWDOWN", "SMSURF", "SMROOT")
    HELD = ("RAINRATE", "LAI")

    def __init__(self, npz):
        self.t0 = np.datetime64(npz["time"][0], "s")
        dt = (np.asarray(npz["time"], "datetime64[s]") - self.t0).astype(np.int64)
        self.step = int(dt[1] - dt[0])
        if not np.all(np.diff(dt) == self.step):
            raise ValueError("forcing times are not equally spaced")
        self.data = {k: np.asarray(npz[k], dtype=F32) for k in self.LINEAR + self.HELD}

    def at(self, when):
        s = int((np.datetime64(when, "s") - self.t0).astype(np.int64))
        i, rem = divmod(s, self.step)
        out = {}
        if rem == 0:
            for k, a in self.data.items():
                out[k] = a[i]
            return out
        idts, idts2 = rem, self.step
        frac = F32(F32(idts2 - idts) / F32(idts2))
        for k in self.LINEAR:
            a = self.data[k]
            out[k] = F32(a[i] * frac + a[i + 1] * (F32(1.0) - frac))
        for k in self.HELD:
            out[k] = self.data[k][i]
        return out


def calc_declin(when, lat, lon):
    """Ports CALC_DECLIN (module_NoahMP_hrldas_driver.F): cosine of the solar zenith angle at time `when` (UTC)."""
    degrad = f32(3.14159265 / 180.0)
    dpd = f32(360.0 / 365.0)
    iday = (when.date() - when.date().replace(month=1, day=1)).days
    gmt = when.hour + when.minute / 60.0 + when.second / 3600.0
    julian = iday + gmt / 24.0
    obecl = 23.5 * degrad
    sinob = math.sin(obecl)
    sxlong = dpd * (julian - 80.0) * degrad if julian >= 80.0 else dpd * (julian + 285.0) * degrad
    declin = math.asin(sinob * math.sin(sxlong))
    tloc = when.hour + when.minute / 60.0 + when.second / 3600.0 + lon / 15.0
    tloc = math.fmod(tloc + 24.0, 24.0)
    hrang = 15.0 * (tloc - 12.0) * degrad
    cosz = math.sin(lat * degrad) * math.sin(declin) + math.cos(lat * degrad) * math.cos(declin) * math.cos(hrang)
    return f32(cosz), julian


# ---------------------------------------------------------------------------------------------------------------------
# point model
# ---------------------------------------------------------------------------------------------------------------------
class PointModel:
    """One HRLDAS point: `io` holds the NoahmpIO fields that persist between steps (TVXY, TGXY, TSLB, ...)."""

    def __init__(self, run_dir):
        self.run_dir = run_dir
        self.nl = tables.read_hrldas_namelist(os.path.join(run_dir, "namelist.hrldas"))
        self.setup = tables.read_setup(self.nl["HRLDAS_SETUP_FILE"])
        if self.setup["MMINLU"].strip() != "MODIFIED_IGBP_MODIS_NOAH":
            raise NotImplementedError("only the MODIS land-use table is ported")
        self.tbl = tables.read_table(os.path.join(run_dir, "NoahmpTable.TBL"), self.setup["IVGTYP"],
                                     self.setup["ISLTYP"])
        self.soil = tables.read_soil_spatial(self.nl["SPATIAL_FILENAME"])
        self.p = build_params(self.nl, self.tbl, self.setup, self.soil)
        self.forcing = HourlyForcing(np.load(os.path.join(os.path.dirname(self.nl["HRLDAS_SETUP_FILE"]),
                                                          "forcing_arrays.npz")))
        nl = self.nl
        self.start = datetime(int(nl["START_YEAR"]), int(nl["START_MONTH"]), int(nl["START_DAY"]),
                              int(nl["START_HOUR"]), int(nl["START_MIN"]))
        self.dt = int(nl["NOAH_TIMESTEP"])
        khour = int(nl["KHOUR"]) if "KHOUR" in nl else int(nl["KDAY"]) * 24
        self.ntime = int(khour * 3600 / self.dt)
        self.init_state()

    # ---------------------------------------------------------------------------------------------------------------
    def init_state(self):
        """Cold start: READINIT_HRLDAS / READVEG_HRLDAS values, then NoahmpInitMain and NoahmpSnowInitMain."""
        S, p, tbl = self.setup, self.p, self.tbl
        n = p.NumSoilLayer
        if not np.allclose(S["DZS"][:n], self.nl["SOIL_THICK_INPUT"], atol=1e-6):
            raise NotImplementedError("setup soil layers differ from the namelist layers (init_interp not ported)")
        io = SimpleNamespace()
        snow = S["SNOW"]
        snowh = S["SNODEP"] if S["SNODEP"] is not None else snow * 0.005
        snow, snowh = max(snow, 0.0), max(snowh, 0.0)
        if snow > 0.0 or snowh > 0.0:
            raise NotImplementedError("initial snow is not ported")
        io.SNOW, io.SNOWH = snow, snowh
        smois = list(S["SMOIS"][:n])
        tslb = list(S["TSLB"][:n])
        bexp, smcmax, psisat = tbl["BEXP_TABLE"], tbl["SMCMAX_TABLE"], tbl["PSISAT_TABLE"]
        smois = [min(x, smcmax) for x in smois]
        sh2o = []
        for k in range(n):
            if bexp > 0.0 and smcmax > 0.0 and psisat > 0.0 and tslb[k] < 273.149:
                fk = (((3.335e5 / (9.81 * (-psisat))) * ((tslb[k] - 273.15) / tslb[k])) ** (-1 / bexp)) * smcmax
                sh2o.append(min(max(fk, 0.02), smois[k]))
            else:
                sh2o.append(smois[k])
        io.SMOIS, io.SH2O, io.TSLB = smois, sh2o, tslb
        tsk = S["TSK"]
        io.TV = io.TG = io.TAH = tsk
        io.CANLIQ = io.CANICE = 0.0
        io.EAH = 2000.0
        io.CM = io.CH = 0.0
        io.FWET = 0.0
        io.SNEQVO = 0.0
        io.TAUSS = 0.0
        lai = S["LAI"]
        if lai > 100 or lai < 0:
            lai = 0.0
        io.LAI = max(lai, 0.05)
        io.XSAI = max(0.1 * io.LAI, 0.05)
        io.ISNOW = 0
        io.ZSNSO = list(p.DepthSoilLayer)
        io.TMN = S["TMN"]
        io.LAT, io.LON = S["XLAT"], S["XLONG"]
        self.io = io

    # ---------------------------------------------------------------------------------------------------------------
    def read_forcing(self, itime):
        """land_driver_exe: READFORC_HRLDAS at OLDDATE = start + itime * dt and the unit conversions of the driver."""
        when = self.start + timedelta(seconds=itime * self.dt)
        f = self.forcing.at(when)
        io = self.io
        io.when = when
        io.T, io.U, io.V, io.P = f["T2D"], f["U2D"], f["V2D"], f["PSFC"]
        io.QV = F32(f["Q2D"] / (F32(1.0) - f["Q2D"]))              # specific humidity -> mixing ratio (driver)
        io.GLW, io.SWDOWN = f["LWDOWN"], f["SWDOWN"]
        io.RAINBL = F32(f["RAINRATE"] * F32(self.dt))              # mm per time step
        io.LAI = f["LAI"]                                          # update_lai (DYNAMIC_VEG_OPTION 8)
        io.SMSURF_PRESC, io.SMROOT_PRESC = f["SMSURF"], f["SMROOT"]
        io.COSZEN, io.JULIAN = calc_declin(when, io.LAT, io.LON)
        return when

    def prescribe_soil_moisture(self):
        """[noah/ s*] module_NoahMP_hrldas_driver.F, 'Prescribe soil moisture': layer 1 <- SMSURF, layers 2..NSOIL <-
        SMROOT; total water replaced, existing ice kept (capped at the new total), liquid = total - ice."""
        io = self.io
        for k in range(self.p.NumSoilLayer):
            sm = float(io.SMSURF_PRESC) if k == 0 else float(io.SMROOT_PRESC)
            if 0.0 <= sm <= 1.0:
                ice = max(0.0, io.SMOIS[k] - io.SH2O[k])
                ice = min(ice, sm)
                io.SMOIS[k] = sm
                io.SH2O[k] = sm - ice

    def column_vars(self):
        """ForcingVarInTransfer, ConfigVarInTransfer, EnergyVarInTransfer, WaterVarInTransfer, BiochemVarInTransfer:
        the per-step column variables built from the persistent NoahmpIO fields."""
        io, p = self.io, self.p
        v = SimpleNamespace()
        # config
        v.NumSnowLayerNeg = io.ISNOW
        v.DepthSnowSoilLayer = list(io.ZSNSO)
        v.CosSolarZenithAngle = float(io.COSZEN)
        # forcing
        v.TemperatureAirRefHeight = float(io.T)
        v.WindEastwardRefHeight = float(io.U)
        v.WindNorthwardRefHeight = float(io.V)
        v.SpecHumidityRefHeight = float(F32(io.QV / (F32(1.0) + io.QV)))
        v.PressureAirRefHeight = float(F32((io.P + io.P) * F32(0.5)))
        v.PressureAirSurface = float(io.P)
        v.RadLwDownRefHeight = float(io.GLW)
        v.RadSwDownRefHeight = float(io.SWDOWN)
        v.TemperatureSoilBottom = io.TMN
        pr = float(F32(io.RAINBL / F32(self.dt)))
        v.PrecipConvRefHeight = 0.0
        v.PrecipNonConvRefHeight = pr
        v.PrecipShConvRefHeight = 0.0
        v.PrecipSnowRefHeight = v.PrecipGraupelRefHeight = v.PrecipHailRefHeight = 0.0
        v.RadSwDirFrac, v.RadSwVisFrac = f32(0.7), f32(0.5)
        # energy state
        v.LeafAreaIndex = float(io.LAI)
        v.StemAreaIndex = io.XSAI
        v.TemperatureGrd = io.TG
        v.TemperatureCanopy = io.TV
        v.SnowAgeNondim = io.TAUSS
        v.PressureVaporCanAir = io.EAH
        v.TemperatureCanopyAir = io.TAH
        v.ExchCoeffShSfc = io.CH
        v.ExchCoeffMomSfc = io.CM
        v.TemperatureSoilSnow = list(io.TSLB)
        v.PressureAtmosCO2 = p.CO2 * v.PressureAirRefHeight
        v.PressureAtmosO2 = p.O2 * v.PressureAirRefHeight
        v.HeatGroundTotAcc = 0.0                                   # ACC_SSOILXY reset (SOIL_UPDATE_STEPS = 1)
        # water state
        v.CanopyLiqWater = io.CANLIQ
        v.CanopyIce = io.CANICE
        v.CanopyWetFrac = io.FWET
        v.SnowWaterEquiv = io.SNOW
        v.SnowWaterEquivPrev = io.SNEQVO
        v.SnowDepth = io.SNOWH
        v.SoilLiqWater = list(io.SH2O)
        v.SoilMoisture = list(io.SMOIS)
        v.SoilMoistureRootPresc = float(io.SMROOT_PRESC) if p.PrescribeSoilMoisture == 1 else -9999.0
        # biochem state
        v.NitrogenConcFoliage = 1.0
        return v

    def store_vars(self, v):
        """The *VarOutTransfer modules: column results back into the persistent NoahmpIO fields."""
        io = self.io
        io.TSLB = list(v.TemperatureSoilSnow)
        io.TV, io.TG = v.TemperatureCanopy, v.TemperatureGrd
        io.EAH, io.TAH = v.PressureVaporCanAir, v.TemperatureCanopyAir
        io.CM, io.CH = v.ExchCoeffMomSfc, v.ExchCoeffShSfc
        io.LAI, io.XSAI = v.LeafAreaIndex, v.StemAreaIndex
        io.TAUSS = v.SnowAgeNondim
        io.SMOIS, io.SH2O = list(v.SoilMoisture), list(v.SoilLiqWater)
        io.CANLIQ, io.CANICE, io.FWET = v.CanopyLiqWater, v.CanopyIce, v.CanopyWetFrac
        io.SNOW, io.SNOWH, io.SNEQVO = v.SnowWaterEquiv, v.SnowDepth, v.SnowWaterEquivPrev
        io.ISNOW = v.NumSnowLayerNeg
        io.ZSNSO = list(v.DepthSnowSoilLayer)

    # ---------------------------------------------------------------------------------------------------------------
    def step(self, itime):
        """land_driver_exe for itime >= 1: forcing at OLDDATE = start + itime*dt (the END of the step), first-step
        guesses, prescription, NoahmpDriverMain. Returns (time stamp, column variables)."""
        when = self.read_forcing(itime)
        io = self.io
        if itime == 1:
            io.EAH = float(F32(io.P * io.QV / (F32(0.622) + io.QV)))
            io.TAH = float(io.T)
            io.CH = io.CM = 0.1
        if self.p.PrescribeSoilMoisture == 1:
            self.prescribe_soil_moisture()
        v = self.column_vars()
        noahmp_main(v, self.p)
        self.store_vars(v)
        return when, v

    def run(self, nsteps=None, record=None):
        """Run steps 1..nsteps (default NTIME). `record(v, io)` -> dict of output values per step.
        Returns (times, {name: array})."""
        nsteps = self.ntime if nsteps is None else nsteps
        self.read_forcing(0)                                       # ITIME = 0: forcing read, no physics
        times, out = [], {}
        for it in range(1, nsteps + 1):
            when, v = self.step(it)
            times.append(when)
            row = (record or ldasout_record)(v, self.io)
            for k, val in row.items():
                out.setdefault(k, []).append(val)
        return np.array(times, dtype="datetime64[s]"), {k: np.array(a, dtype=np.float64) for k, a in out.items()}


def ldasout_record(v, io):
    """Values as written to LDASOUT (EnergyVarOutTransfer, WaterVarOutTransfer, driver forcing fields)."""
    sh2m = lambda q: q / (1.0 - q)
    return {
        "FSA": v.RadSwAbsSfc, "FIRA": v.RadLwNetSfc, "SAV": v.RadSwAbsVeg, "SAG": v.RadSwAbsGrd,
        "ALBEDO": v.AlbedoSfc, "EMISS": v.EmissivitySfc, "COSZ": float(io.COSZEN), "FVEG": v.VegFrac,
        "LAI": v.LeafAreaIndex, "SAI": v.StemAreaIndex, "TRAD": v.TemperatureRadSfc, "TV": v.TemperatureCanopy,
        "TG": v.TemperatureGrd, "TGV": v.TemperatureGrdVeg, "TGB": v.TemperatureGrdBare,
        "TAH": v.TemperatureCanopyAir, "EAH": v.PressureVaporCanAir, "T2MV": v.TemperatureAir2mVeg,
        "T2MB": v.TemperatureAir2mBare, "Q2MV": sh2m(v.SpecHumidity2mVeg), "Q2MB": sh2m(v.SpecHumidity2mBare),
        "LH": v.HeatLatentGrd + v.HeatLatentCanopy + v.HeatLatentTransp, "HFX": v.HeatSensibleSfc,
        "GRDFLX": v.HeatGroundTot, "TR": v.HeatLatentCanTransp, "EVC": v.HeatLatentCanEvap,
        "IRC": v.RadLwNetCanopy, "SHC": v.HeatSensibleCanopy, "IRG": v.RadLwNetVegGrd, "SHG": v.HeatSensibleVegGrd,
        "EVG": v.HeatLatentVegGrd, "GHV": v.HeatGroundVegGrd, "IRB": v.RadLwNetBareGrd,
        "SHB": v.HeatSensibleBareGrd, "EVB": v.HeatLatentBareGrd, "GHB": v.HeatGroundBareGrd,
        "ECAN": v.EvapCanopyNet, "ETRAN": v.Transpiration, "EDIR": v.EvapGroundNet, "CANLIQ": v.CanopyLiqWater,
        "FWET": v.CanopyWetFrac, "CM": v.ExchCoeffMomSfc, "CH": v.ExchCoeffShSfc, "CHV": v.ExchCoeffShAbvCan,
        "CHB": v.ExchCoeffShBare, "CHLEAF": v.ExchCoeffShLeaf, "CHUC": v.ExchCoeffShUndCan,
        "CHV2": v.ExchCoeffSh2mVeg, "CHB2": v.ExchCoeffSh2mBare, "BTRAN": v.SoilTranspFacAcc,
        "SWFORC": float(io.SWDOWN), "LWFORC": float(io.GLW), "RAINRATE": float(io.RAINBL),
        "APAR": v.RadPhotoActAbsCan, "PSN": v.PhotosynTotal, "CANHS": v.HeatCanStorageChg,
        "SOIL_T1": v.TemperatureSoilSnow[0], "SOIL_T2": v.TemperatureSoilSnow[1],
        "SOIL_T3": v.TemperatureSoilSnow[2], "SOIL_T4": v.TemperatureSoilSnow[3],
        "SOIL_M1": v.SoilMoisture[0], "SOIL_M2": v.SoilMoisture[1],
        "SOIL_W1": v.SoilLiqWater[0], "SOIL_W2": v.SoilLiqWater[1],
    }
