"""Atmospheric forcing processing, general initialization and phenology (one column, one time step).

Every function takes `v` (per-step variables, Fortran names) and `p` (parameters / options) and updates `v` in place.
"""
import math

from .constants import RAIR, CPAIR, TFRZ


def process_atmos_forcing(v, p):
    """Ports noahmp/src/AtmosForcingMod.F90 (ProcessAtmosForcing). Rain/snow partition: option 1 (Jordan 1991) only."""
    T, P, q = v.TemperatureAirRefHeight, v.PressureAirRefHeight, v.SpecHumidityRefHeight
    v.TemperaturePotRefHeight = T * (P / P) ** (RAIR / CPAIR)
    v.PressureVaporRefHeight = q * P / (0.622 + 0.378 * q)
    v.DensityAirRefHeight = (P - 0.378 * v.PressureVaporRefHeight) / (RAIR * T)

    # downward solar radiation (direct / diffuse, visible / NIR)
    rad_dir, rad_vis = 0.7, 0.5
    if 0.0 <= v.RadSwDirFrac <= 1.0:
        rad_dir = v.RadSwDirFrac
    else:
        v.RadSwDirFrac = rad_dir
    if 0.0 <= v.RadSwVisFrac <= 1.0:
        rad_vis = v.RadSwVisFrac
    else:
        v.RadSwVisFrac = rad_vis
    if v.CosSolarZenithAngle <= 0.0:
        v.RadSwDownRefHeight = 0.0
    sw = v.RadSwDownRefHeight
    v.RadSwDownDir = [sw * rad_dir * rad_vis, sw * rad_dir * (1.0 - rad_vis)]
    v.RadSwDownDif = [sw * (1.0 - rad_dir) * rad_vis, sw * (1.0 - rad_dir) * (1.0 - rad_vis)]

    # precipitation
    v.PrecipTotRefHeight = v.PrecipConvRefHeight + v.PrecipNonConvRefHeight + v.PrecipShConvRefHeight
    if p.OptRainSnowPartition == 4:
        raise NotImplementedError("PCP_PARTITION_OPTION 4 is not ported")
    v.PrecipConvTotRefHeight = 0.10 * v.PrecipTotRefHeight
    v.PrecipLargeSclRefHeight = 0.90 * v.PrecipTotRefHeight
    v.PrecipAreaFrac = 0.0
    if (v.PrecipConvTotRefHeight + v.PrecipLargeSclRefHeight) > 0.0:
        v.PrecipAreaFrac = (v.PrecipConvTotRefHeight + v.PrecipLargeSclRefHeight) / \
                           (10.0 * v.PrecipConvTotRefHeight + v.PrecipLargeSclRefHeight)
    if p.OptRainSnowPartition != 1:
        raise NotImplementedError("only PCP_PARTITION_OPTION 1 (Jordan 1991) is ported")
    if T > (TFRZ + 2.5):
        v.FrozenPrecipFrac = 0.0
    elif T <= (TFRZ + 0.5):
        v.FrozenPrecipFrac = 1.0
    elif T <= (TFRZ + 2.0):
        v.FrozenPrecipFrac = 1.0 - (-54.632 + 0.2 * T)
    else:
        v.FrozenPrecipFrac = 0.6
    v.SnowfallDensity = min(p.SnowfallDensityMax, 67.92 + 51.25 * math.exp((T - TFRZ) / 2.59))
    v.RainfallRefHeight = v.PrecipTotRefHeight * (1.0 - v.FrozenPrecipFrac)
    v.SnowfallRefHeight = v.PrecipTotRefHeight * v.FrozenPrecipFrac

    # wind speed at reference height for turbulence calculation
    v.WindSpdRefHeight = max(math.sqrt(v.WindEastwardRefHeight ** 2 + v.WindNorthwardRefHeight ** 2), 1.0)


def general_init(v, p):
    """Ports noahmp/src/GeneralInitMod.F90 (GeneralInit): layer thickness, root-zone temperature, land-type flags.
    No snow layers (NumSnowLayerNeg = 0), so the layers are the soil layers."""
    z = v.DepthSnowSoilLayer
    v.ThicknessSnowSoilLayer = [-z[0]] + [z[k - 1] - z[k] for k in range(1, p.NumSoilLayer)]
    nr = p.NumSoilLayerRoot
    v.TemperatureRootZone = 0.0
    for k in range(nr):
        v.TemperatureRootZone += v.TemperatureSoilSnow[k] * v.ThicknessSnowSoilLayer[k] / (-p.DepthSoilLayer[nr - 1])
    # MODIFIED_IGBP_MODIS_NOAH
    v.FlagCropland = p.VegType in (12, 14)
    v.FlagWetland = p.VegType == 11


def phenology_main(v, p):
    """Ports noahmp/src/PhenologyMainMod.F90 (PhenologyMain) for DYNAMIC_VEG_OPTION 8 (LAI from input, vegetation
    fraction from LAI + SAI), no crop."""
    if p.OptDynamicVeg != 8:
        raise NotImplementedError("only DYNAMIC_VEG_OPTION 8 is ported")
    lai = v.LeafAreaIndex
    sai = max(0.05, 0.1 * lai)
    if lai < 0.05:
        sai = 0.0
    if sai < 0.05:
        sai = 0.0
    if (lai < 0.05) or (sai == 0.0):
        lai = 0.0
    if p.VegType in (p.IndexWaterPoint, p.IndexBarrenPoint, p.IndexIcePoint) or p.FlagUrban:
        lai, sai = 0.0, 0.0
    v.LeafAreaIndex, v.StemAreaIndex = lai, sai

    # vegetation fraction buried by snow
    hvt, hvb, snowh = p.HeightCanopyTop, p.HeightCanopyBot, v.SnowDepth
    thick_bury = min(max(snowh - hvb, 0.0), (hvt - hvb))
    v.CanopyFracSnowBury = thick_bury / max(1.0e-06, (hvt - hvb))
    if (hvt > 0.0) and (hvt <= 1.0):
        snow_veg_bury = hvt * math.exp(-min(snowh, 10.0) / 0.2)
        v.CanopyFracSnowBury = min(snowh, snow_veg_bury) / snow_veg_bury

    v.LeafAreaIndEff = lai * (1.0 - v.CanopyFracSnowBury)
    v.StemAreaIndEff = sai * (1.0 - v.CanopyFracSnowBury)
    if v.StemAreaIndEff < 0.05:
        v.StemAreaIndEff = 0.0
    if (v.LeafAreaIndEff < 0.05) or (v.StemAreaIndEff == 0.0):
        v.LeafAreaIndEff = 0.0

    v.IndexGrowSeason = 1.0 if v.TemperatureCanopy > p.TemperatureMinPhotosyn else 0.0

    v.VegFrac = 1.0 - math.exp(-0.52 * (lai + sai))
    if v.VegFrac <= 0.05:
        v.VegFrac = 0.05
    if p.FlagUrban or (p.VegType == p.IndexBarrenPoint):
        v.VegFrac = 0.0
    if (v.LeafAreaIndEff + v.StemAreaIndEff) == 0.0:
        v.VegFrac = 0.0
