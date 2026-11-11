"""Ground properties used by the energy balance: roughness, soil thermal properties, soil water transpiration factor
(BTRAN), soil evaporation resistance and psychrometric variables (one column, one time step, no snow layers)."""
import math

from .constants import (CWAT, CICE, CPAIR, TKWAT, TKICE, TKQTZ, TKO, TFRZ, HVAP, HSUB, GRAV, RW)


def ground_roughness_property(v, p, flag_veg_sfc):
    """Ports noahmp/src/GroundRoughnessPropertyMod.F90 (GroundRoughnessProperty), soil point, not urban."""
    v.RoughLenMomGrd = p.RoughLenMomSoil * (1.0 - v.SnowCoverFrac) + v.SnowCoverFrac * p.RoughLenMomSnow
    v.ZeroPlaneDispGrd = v.SnowDepth
    if flag_veg_sfc:
        v.RoughLenMomSfc = p.RoughLenMomVeg
        v.ZeroPlaneDispSfc = 0.65 * p.HeightCanopyTop
        if v.SnowDepth > v.ZeroPlaneDispSfc:
            v.ZeroPlaneDispSfc = v.SnowDepth
    else:
        v.RoughLenMomSfc = v.RoughLenMomGrd
        v.ZeroPlaneDispSfc = v.ZeroPlaneDispGrd
    v.RefHeightAboveGrd = v.ZeroPlaneDispSfc + p.RefHeightAboveSfc
    if v.ZeroPlaneDispGrd >= v.RefHeightAboveGrd:
        v.RefHeightAboveGrd = v.ZeroPlaneDispGrd + p.RefHeightAboveSfc


def soil_thermal_property(v, p):
    """Ports noahmp/src/SoilThermalPropertyMod.F90 (SoilThermalProperty): heat capacity and thermal conductivity
    (Johansen / Peters-Lidard, 'fine' Kersten number)."""
    n = p.NumSoilLayer
    v.HeatCapacVolSoil = [0.0] * n
    v.ThermConductSoil = [0.0] * n
    for k in range(n):
        smc, sh2o, smcmax = v.SoilMoisture[k], v.SoilLiqWater[k], p.SoilMoistureSat[k]
        ice = smc - sh2o
        v.HeatCapacVolSoil[k] = sh2o * CWAT + (1.0 - smcmax) * p.SoilHeatCapacity + (smcmax - smc) * CPAIR + ice * CICE
        sat_ratio = smc / smcmax
        qz = p.SoilQuartzFrac[k]
        tk_solid = (TKQTZ ** qz) * (TKO ** (1.0 - qz))
        wat_frac = 1.0
        if smc > 0.0:
            wat_frac = sh2o / smc
        wat_frac_sat = wat_frac * smcmax
        tk_sat = tk_solid ** (1.0 - smcmax) * TKICE ** (smcmax - wat_frac_sat) * TKWAT ** wat_frac_sat
        gam = (1.0 - smcmax) * 2700.0
        tk_dry = (0.135 * gam + 64.7) / (2700.0 - 0.947 * gam)
        if (sh2o + 0.0005) < smc:          # frozen
            kersten = sat_ratio
        else:
            kersten = math.log10(sat_ratio) + 1.0 if sat_ratio > 0.1 else 0.0
        v.ThermConductSoil[k] = kersten * (tk_sat - tk_dry) + tk_dry


def ground_thermal_property(v, p):
    """Ports noahmp/src/GroundThermalPropertyMod.F90 (GroundThermalProperty) without snow layers
    (SnowThermalProperty has no layer to act on), soil point, not urban."""
    if v.NumSnowLayerNeg != 0:
        raise NotImplementedError("snow layers are not ported")
    soil_thermal_property(v, p)
    n = p.NumSoilLayer
    dz = v.ThicknessSnowSoilLayer
    v.ThermConductSoilSnow = list(v.ThermConductSoil)
    v.HeatCapacSoilSnow = list(v.HeatCapacVolSoil)
    v.PhaseChgFacSoilSnow = [p.MainTimeStep / (v.HeatCapacSoilSnow[k] * dz[k]) for k in range(n)]
    # snow/soil interface without snow layer
    v.ThermConductSoilSnow[0] = (v.ThermConductSoilSnow[0] * dz[0] + 0.35 * v.SnowDepth) / (v.SnowDepth + dz[0])


def soil_water_transpiration(v, p):
    """Ports noahmp/src/SoilWaterTranspirationMod.F90 (SoilWaterTranspiration), BTR_OPTION 4 ([noah/ s*] Laio beta on
    one root-zone soil moisture: the prescribed SMROOT when in [0, 1], else the thickness-weighted root-layer liquid
    water)."""
    if p.OptSoilWaterTranspiration != 4:
        raise NotImplementedError("only BTR_OPTION 4 is ported")
    nr = p.NumSoilLayerRoot
    dz = v.ThicknessSnowSoilLayer
    zr = -p.DepthSoilLayer[nr - 1]
    if 0.0 <= v.SoilMoistureRootPresc <= 1.0:
        theta = v.SoilMoistureRootPresc
    else:
        theta = 0.0
        for k in range(nr):
            theta += v.SoilLiqWater[k] * dz[k] / zr
    btr = min(1.0, max(0.0, (theta - p.SoilMoistureLaioWilt) / (p.SoilMoistureLaioStar - p.SoilMoistureLaioWilt)))
    if getattr(p, "BtranOverride", None) is not None:      # sensitivity hook (not in Fortran): fixed beta, ignores soil moisture
        btr = float(p.BtranOverride)
    v.SoilTranspFacAcc = max(1.0e-6, btr)
    v.SoilTranspFac = [dz[k] / zr for k in range(nr)]


def resistance_ground_evaporation(v, p):
    """Ports noahmp/src/ResistanceGroundEvaporationMod.F90 (ResistanceGroundEvaporation), soil point,
    SURFACE_RESISTANCE_OPTION 1 or 4 (Sakaguchi & Zeng 2009; option 4 adds the snow-cover weighting), with the
    [noah/ s*] nuisance scale."""
    if p.OptGroundResistanceEvap not in (1, 4):
        raise NotImplementedError("only SURFACE_RESISTANCE_OPTION 1 / 4 are ported")
    sh2o1, smcmax1 = v.SoilLiqWater[0], p.SoilMoistureSat[0]
    dry = (-p.DepthSoilLayer[0]) * (math.exp((1.0 - min(1.0, sh2o1 / smcmax1)) ** p.ResistanceSoilExp) - 1.0) / \
        (2.71828 - 1.0)
    vap_red = 2.2e-5 * smcmax1 * smcmax1 * (1.0 - p.SoilMoistureWilt[0] / smcmax1) ** (2.0 + 3.0 / p.SoilExpCoeffB[0])
    rs = dry / vap_red
    rs = rs * p.ResistanceGrdEvapScale
    scf = v.SnowCoverFrac
    if p.OptGroundResistanceEvap == 4:
        rs = 1.0 / (scf * (1.0 / p.ResistanceSnowSfc) + (1.0 - scf) * (1.0 / max(rs, 0.001)))
    if (sh2o1 < 0.01) and (v.SnowDepth == 0.0):
        rs = 1.0e6
    v.ResistanceGrdEvap = rs
    psi = -p.SoilMatPotentialSat[0] * (max(0.01, sh2o1) / smcmax1) ** (-p.SoilExpCoeffB[0])
    v.RelHumidityGrd = scf + (1.0 - scf) * math.exp(psi * GRAV / (RW * v.TemperatureGrd))


def psychrometric_variable(v, p):
    """Ports noahmp/src/PsychrometricVariableMod.F90 (PsychrometricVariable)."""
    P = v.PressureAirRefHeight
    if v.TemperatureCanopy > TFRZ:
        v.LatHeatVapCanopy, v.FlagFrozenCanopy = HVAP, False
    else:
        v.LatHeatVapCanopy, v.FlagFrozenCanopy = HSUB, True
    v.PsychConstCanopy = CPAIR * P / (0.622 * v.LatHeatVapCanopy)
    if v.TemperatureGrd > TFRZ:
        v.LatHeatVapGrd, v.FlagFrozenGround = HVAP, False
    else:
        v.LatHeatVapGrd, v.FlagFrozenGround = HSUB, True
    v.PsychConstGrd = CPAIR * P / (0.622 * v.LatHeatVapGrd)
