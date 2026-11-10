"""Canopy interception, precipitation heat advection and canopy hydrology (one column, one time step)."""
import math

from .constants import CWAT, CICE, HVAP, HSUB, HFUS, TFRZ, DENH2O, DENICE


def canopy_water_intercept(v, p):
    """Ports noahmp/src/CanopyWaterInterceptMod.F90 (CanopyWaterIntercept)."""
    dt = p.MainTimeStep
    vai = v.LeafAreaIndEff + v.StemAreaIndEff
    fveg = v.VegFrac
    v.InterceptCanopyRain = v.DripCanopyRain = v.ThroughfallRain = 0.0
    v.InterceptCanopySnow = v.DripCanopySnow = v.ThroughfallSnow = 0.0
    v.RainfallGround = v.SnowfallGround = v.SnowDepthIncr = 0.0
    snow_drip = 0.0

    # canopy liquid water
    v.CanopyLiqWaterMax = fveg * p.CanopyLiqHoldCap * vai
    rain = v.RainfallRefHeight
    if vai > 0.0:
        intr = fveg * rain * v.PrecipAreaFrac
        intr = min(intr, (v.CanopyLiqWaterMax - v.CanopyLiqWater) / dt *
                   (1.0 - math.exp(-rain * dt / v.CanopyLiqWaterMax)))
        intr = max(intr, 0.0)
        v.InterceptCanopyRain = intr
        v.DripCanopyRain = fveg * rain - intr
        v.ThroughfallRain = (1.0 - fveg) * rain
        v.CanopyLiqWater = max(0.0, v.CanopyLiqWater + intr * dt)
    else:
        v.ThroughfallRain = rain
        if v.CanopyLiqWater > 0.0:
            v.DripCanopyRain += v.CanopyLiqWater / dt
            v.CanopyLiqWater = 0.0

    # canopy ice
    v.CanopyIceMax = fveg * 6.6 * (0.27 + 46.0 / v.SnowfallDensity) * vai
    snow = v.SnowfallRefHeight
    if vai > 0.0:
        ints = fveg * snow * v.PrecipAreaFrac
        ints = min(ints, (v.CanopyIceMax - v.CanopyIce) / dt * (1.0 - math.exp(-snow * dt / v.CanopyIceMax)))
        ints = max(ints, 0.0)
        drip_t = max(0.0, (v.TemperatureCanopy - 270.15) / 1.87e5)
        drip_w = math.sqrt(v.WindEastwardRefHeight ** 2 + v.WindNorthwardRefHeight ** 2) / 1.56e5
        snow_drip = max(0.0, v.CanopyIce) * (drip_w + drip_t)
        snow_drip = min(v.CanopyIce / dt + ints, snow_drip)
        v.InterceptCanopySnow = ints
        v.DripCanopySnow = (fveg * snow - ints) + snow_drip
        v.ThroughfallSnow = (1.0 - fveg) * snow
        v.CanopyIce = max(0.0, v.CanopyIce + (ints - snow_drip) * dt)
    else:
        v.ThroughfallSnow = snow
        if v.CanopyIce > 0.0:
            v.DripCanopySnow += v.CanopyIce / dt
            v.CanopyIce = 0.0

    # wetted fraction of canopy
    if v.CanopyIce > 0.0:
        fwet = max(0.0, v.CanopyIce) / max(v.CanopyIceMax, 1.0e-06)
    else:
        fwet = max(0.0, v.CanopyLiqWater) / max(v.CanopyLiqWaterMax, 1.0e-06)
    v.CanopyWetFrac = min(fwet, 1.0) ** 0.667
    v.CanopyTotalWater = v.CanopyLiqWater + v.CanopyIce

    # rain or snow on the ground
    v.RainfallGround = v.DripCanopyRain + v.ThroughfallRain
    v.SnowfallGround = v.DripCanopySnow + v.ThroughfallSnow
    v.SnowDepthIncr = v.SnowfallGround / v.SnowfallDensity


def precipitation_heat_advect(v, p):
    """Ports noahmp/src/PrecipitationHeatAdvectMod.F90 (PrecipitationHeatAdvect)."""
    ta, tv, tg, fveg = v.TemperatureAirRefHeight, v.TemperatureCanopy, v.TemperatureGrd, v.VegFrac
    cw, ci = CWAT / 1000.0, CICE / 1000.0
    air_can = fveg * v.RainfallRefHeight * cw * (ta - tv)
    can_grd = v.DripCanopyRain * cw * (tv - tg)
    air_grd = v.ThroughfallRain * cw * (ta - tg)
    air_can += fveg * v.SnowfallRefHeight * ci * (ta - tv)
    can_grd += v.DripCanopySnow * ci * (tv - tg)
    air_grd += v.ThroughfallSnow * ci * (ta - tg)

    pah_v = air_can - can_grd
    pah_g = can_grd
    pah_b = air_grd
    if 0.0 < fveg < 1.0:
        pah_g = pah_g / fveg
        pah_b = pah_b / (1.0 - fveg)
    elif fveg <= 0.0:
        pah_b = pah_g + pah_b
        pah_g = 0.0
        pah_v = 0.0
    elif fveg >= 1.0:
        pah_b = 0.0
    v.HeatPrecipAdvCanopy = min(max(pah_v, -20.0), 20.0)
    v.HeatPrecipAdvVegGrd = min(max(pah_g, -20.0), 20.0)
    v.HeatPrecipAdvBareGrd = min(max(pah_b, -20.0), 20.0)


def canopy_hydrology(v, p):
    """Ports noahmp/src/CanopyHydrologyMod.F90 (CanopyHydrology): canopy evaporation / dew / sublimation / frost and
    canopy phase change, after the energy solution."""
    dt = p.MainTimeStep
    vai = v.LeafAreaIndEff + v.StemAreaIndEff
    v.EvapCanopyNet = v.Transpiration = v.EvapCanopyLiq = v.DewCanopyLiq = 0.0
    v.FrostCanopyIce = v.SublimCanopyIce = v.MeltCanopyIce = v.FreezeCanopyLiq = 0.0

    v.CanopyLiqWaterMax = v.VegFrac * p.CanopyLiqHoldCap * vai
    if not v.FlagFrozenCanopy:
        v.Transpiration = max(v.HeatLatentTransp / HVAP, 0.0)
        v.EvapCanopyLiq = max(v.HeatLatentCanopy / HVAP, 0.0)
        v.DewCanopyLiq = abs(min(v.HeatLatentCanopy / HVAP, 0.0))
    else:
        v.Transpiration = max(v.HeatLatentTransp / HSUB, 0.0)
        v.SublimCanopyIce = max(v.HeatLatentCanopy / HSUB, 0.0)
        v.FrostCanopyIce = abs(min(v.HeatLatentCanopy / HSUB, 0.0))

    v.EvapCanopyLiq = min(v.CanopyLiqWater / dt, v.EvapCanopyLiq)
    v.CanopyLiqWater = max(0.0, v.CanopyLiqWater + (v.DewCanopyLiq - v.EvapCanopyLiq) * dt)
    if v.CanopyLiqWater <= 1.0e-06:
        v.CanopyLiqWater = 0.0

    v.CanopyIceMax = v.VegFrac * 6.6 * (0.27 + 46.0 / v.SnowfallDensity) * vai
    v.SublimCanopyIce = min(v.CanopyIce / dt, v.SublimCanopyIce)
    v.CanopyIce = max(0.0, v.CanopyIce + (v.FrostCanopyIce - v.SublimCanopyIce) * dt)
    if v.CanopyIce <= 1.0e-6:
        v.CanopyIce = 0.0

    if (v.CanopyIce > 0.0) and (v.CanopyIce >= v.CanopyLiqWater):
        fwet = max(0.0, v.CanopyIce) / max(v.CanopyIceMax, 1.0e-06)
    else:
        fwet = max(0.0, v.CanopyLiqWater) / max(v.CanopyLiqWaterMax, 1.0e-06)
    v.CanopyWetFrac = min(fwet, 1.0) ** 0.667
    v.CanopyTotalWater = v.CanopyLiqWater + v.CanopyIce

    # canopy ice melting
    if (v.CanopyIce > 1.0e-6) and (v.TemperatureCanopy > TFRZ):
        v.MeltCanopyIce = min(v.CanopyIce / dt, (v.TemperatureCanopy - TFRZ) * CICE * v.CanopyIce / DENICE /
                              (dt * HFUS))
        v.CanopyIce = max(0.0, v.CanopyIce - v.MeltCanopyIce * dt)
        v.CanopyLiqWater = max(0.0, v.CanopyTotalWater - v.CanopyIce)
        v.TemperatureCanopy = v.CanopyWetFrac * TFRZ + (1.0 - v.CanopyWetFrac) * v.TemperatureCanopy
    # canopy water refreezing
    if (v.CanopyLiqWater > 1.0e-6) and (v.TemperatureCanopy < TFRZ):
        v.FreezeCanopyLiq = min(v.CanopyLiqWater / dt, (TFRZ - v.TemperatureCanopy) * CWAT * v.CanopyLiqWater /
                                DENH2O / (dt * HFUS))
        v.CanopyLiqWater = max(0.0, v.CanopyLiqWater - v.FreezeCanopyLiq * dt)
        v.CanopyIce = max(0.0, v.CanopyTotalWater - v.CanopyLiqWater)
        v.TemperatureCanopy = v.CanopyWetFrac * TFRZ + (1.0 - v.CanopyWetFrac) * v.TemperatureCanopy

    v.CanopyTotalWater = v.CanopyLiqWater + v.CanopyIce
    v.EvapCanopyNet = v.EvapCanopyLiq + v.SublimCanopyIce - v.DewCanopyLiq - v.FrostCanopyIce
