"""Column model for one time step: NoahmpMain -> (forcing, phenology, interception, heat advection, EnergyMain,
WaterMain subset)."""
from .constants import SB, UNDEFINED, TFRZ
from . import atmos, canopy_water, radiation, ground, turbulence, soil_thermal


def energy_main(v, p):
    """Ports noahmp/src/EnergyMainMod.F90 (EnergyMain) for a vegetated soil point without snow layers."""
    for name in ("WindStressEwVeg", "WindStressNsVeg", "RadLwNetCanopy", "HeatSensibleCanopy", "RadLwNetVegGrd",
                 "HeatSensibleVegGrd", "HeatLatentVegGrd", "HeatLatentCanEvap", "HeatLatentCanTransp",
                 "HeatGroundVegGrd", "PhotosynLeafSunlit", "PhotosynLeafShade", "TemperatureAir2mVeg",
                 "SpecHumidity2mVeg", "ExchCoeffShAbvCan", "ExchCoeffShLeaf", "ExchCoeffShUndCan", "ExchCoeffSh2mVeg",
                 "HeatPrecipAdvSfc", "HeatCanStorageChg"):
        setattr(v, name, 0.0)

    v.VegAreaIndEff = v.LeafAreaIndEff + v.StemAreaIndEff
    flag_veg_sfc = v.VegAreaIndEff > 0.0

    if p.OptSnowCoverGround != 1:
        raise NotImplementedError("only SNOW_COVER_OPTION 1 is ported")
    radiation.snow_cover_ground_niu07(v, p)
    ground.ground_roughness_property(v, p, flag_veg_sfc)
    ground.ground_thermal_property(v, p)
    radiation.surface_albedo(v, p)
    radiation.surface_radiation(v, p)
    radiation.surface_emissivity(v, p)
    ground.soil_water_transpiration(v, p)
    ground.resistance_ground_evaporation(v, p)
    ground.psychrometric_variable(v, p)

    veg = flag_veg_sfc and (v.VegFrac > 0)
    if veg:
        v.TemperatureGrdVeg = v.TemperatureGrd
        v.ExchCoeffMomAbvCan = v.ExchCoeffMomSfc
        v.ExchCoeffShAbvCan = v.ExchCoeffShSfc
        turbulence.surface_energy_flux_vegetated(v, p)
    v.TemperatureGrdBare = v.TemperatureGrd
    v.ExchCoeffMomBare = v.ExchCoeffMomSfc
    v.ExchCoeffShBare = v.ExchCoeffShSfc
    turbulence.surface_energy_flux_bare_ground(v, p)

    fv = v.VegFrac
    if veg:
        v.WindStressEwSfc = fv * v.WindStressEwVeg + (1.0 - fv) * v.WindStressEwBare
        v.WindStressNsSfc = fv * v.WindStressNsVeg + (1.0 - fv) * v.WindStressNsBare
        v.RadLwNetSfc = fv * v.RadLwNetVegGrd + (1.0 - fv) * v.RadLwNetBareGrd + v.RadLwNetCanopy
        v.HeatSensibleSfc = fv * v.HeatSensibleVegGrd + (1.0 - fv) * v.HeatSensibleBareGrd + v.HeatSensibleCanopy
        v.HeatLatentGrd = fv * v.HeatLatentVegGrd + (1.0 - fv) * v.HeatLatentBareGrd
        v.HeatGroundTot = fv * v.HeatGroundVegGrd + (1.0 - fv) * v.HeatGroundBareGrd
        v.HeatLatentCanopy = v.HeatLatentCanEvap
        v.HeatLatentTransp = v.HeatLatentCanTransp
        v.HeatPrecipAdvSfc = fv * v.HeatPrecipAdvVegGrd + (1.0 - fv) * v.HeatPrecipAdvBareGrd + v.HeatPrecipAdvCanopy
        v.TemperatureGrd = fv * v.TemperatureGrdVeg + (1.0 - fv) * v.TemperatureGrdBare
        v.TemperatureAir2m = fv * v.TemperatureAir2mVeg + (1.0 - fv) * v.TemperatureAir2mBare
        v.TemperatureSfc = fv * v.TemperatureCanopy + (1.0 - fv) * v.TemperatureGrdBare
        v.ExchCoeffMomSfc = fv * v.ExchCoeffMomAbvCan + (1.0 - fv) * v.ExchCoeffMomBare
        v.ExchCoeffShSfc = fv * v.ExchCoeffShAbvCan + (1.0 - fv) * v.ExchCoeffShBare
        v.SpecHumidity2m = fv * v.SpecHumidity2mVeg + (1.0 - fv) * v.SpecHumidity2mBare
        e = v.PressureVaporCanAir
        v.SpecHumiditySfcMean = fv * (e * 0.622 / (v.PressureAirRefHeight - 0.378 * e)) + \
            (1.0 - fv) * v.SpecHumiditySfc
        v.RoughLenMomSfcToAtm = v.RoughLenMomSfc
    else:
        v.WindStressEwSfc = v.WindStressEwBare
        v.WindStressNsSfc = v.WindStressNsBare
        v.RadLwNetSfc = v.RadLwNetBareGrd
        v.HeatSensibleSfc = v.HeatSensibleBareGrd
        v.HeatLatentGrd = v.HeatLatentBareGrd
        v.HeatGroundTot = v.HeatGroundBareGrd
        v.TemperatureGrd = v.TemperatureGrdBare
        v.TemperatureAir2m = v.TemperatureAir2mBare
        v.HeatLatentCanopy = 0.0
        v.HeatLatentTransp = 0.0
        v.HeatPrecipAdvSfc = v.HeatPrecipAdvBareGrd
        v.TemperatureSfc = v.TemperatureGrd
        v.ExchCoeffMomSfc = v.ExchCoeffMomBare
        v.ExchCoeffShSfc = v.ExchCoeffShBare
        v.SpecHumiditySfcMean = v.SpecHumiditySfc
        v.SpecHumidity2m = v.SpecHumidity2mBare
        v.ResistanceStomataSunlit = 0.0
        v.ResistanceStomataShade = 0.0
        v.TemperatureGrdVeg = v.TemperatureGrdBare
        v.ExchCoeffShAbvCan = v.ExchCoeffShBare
        v.RoughLenMomSfcToAtm = v.RoughLenMomGrd

    # emitted longwave and radiative temperature
    v.RadLwEmitSfc = v.RadLwDownRefHeight + v.RadLwNetSfc
    if v.RadLwEmitSfc <= 0.0:
        raise RuntimeError("emitted longwave <= 0")
    es = v.EmissivitySfc
    v.TemperatureRadSfc = ((v.RadLwEmitSfc - (1.0 - es) * v.RadLwDownRefHeight) / (es * SB)) ** 0.25

    v.RadPhotoActAbsCan = v.RadPhotoActAbsSunlit * v.LeafAreaIndSunlit + v.RadPhotoActAbsShade * v.LeafAreaIndShade
    v.PhotosynTotal = v.PhotosynLeafSunlit * v.LeafAreaIndSunlit + v.PhotosynLeafShade * v.LeafAreaIndShade

    # soil temperature at the soil time step (= main time step here: SOIL_TIMESTEP 0)
    v.HeatFromSoilBot = 0.0
    v.HeatGroundTotAcc = v.HeatGroundTotAcc + v.HeatGroundTot
    if p.FlagSoilProcess:
        v.HeatGroundTotMean = v.HeatGroundTotAcc / p.NumSoilTimeStep
        soil_thermal.soil_snow_temperature_main(v, p)

    if p.OptSnowSoilTempTime == 2:
        raise NotImplementedError("TEMP_TIME_SCHEME_OPTION 2 is not ported")

    soil_thermal.soil_snow_water_phase_change(v, p)

    if v.RadSwDownRefHeight > 0.0:
        v.AlbedoSfc = v.RadSwReflSfc / v.RadSwDownRefHeight
    else:
        v.AlbedoSfc = UNDEFINED


def water_main_subset(v, p):
    """Ports the parts of noahmp/src/WaterMainMod.F90 (WaterMain) that feed the next step's energy balance:
    ground evaporation / dew partition, CanopyHydrology, the snow checks of SnowWaterMain (no snow allowed) and the
    frozen-ground surface vapour exchange on the top-layer ice. The soil water solver (SoilWaterMain), runoff and
    groundwater are NOT ported: total soil water is overwritten by the prescription before the next step, and they
    preserve the soil ice, which is the only part of the soil water the prescription keeps."""
    dt = p.MainTimeStep
    n = p.NumSoilLayer
    v.SoilIce = [max(0.0, v.SoilMoisture[k] - v.SoilLiqWater[k]) for k in range(n)]
    v.SnowWaterEquivPrev = v.SnowWaterEquiv
    v.VaporizeGrd = max(v.HeatLatentGrd / v.LatHeatVapGrd, 0.0)
    v.CondenseVapGrd = abs(min(v.HeatLatentGrd / v.LatHeatVapGrd, 0.0))
    v.EvapGroundNet = v.VaporizeGrd - v.CondenseVapGrd

    canopy_water.canopy_hydrology(v, p)

    v.SublimSnowSfcIce = 0.0
    if v.SnowWaterEquiv > 0.0:
        raise NotImplementedError("snow is not ported")
    v.EvapSoilSfcLiq = v.VaporizeGrd - v.SublimSnowSfcIce
    v.FrostSnowSfcIce = 0.0
    v.DewSoilSfcLiq = v.CondenseVapGrd - v.FrostSnowSfcIce

    # SnowWaterMain: with no snow and no snowfall reaching the ground it changes nothing
    if v.SnowfallGround > 0.0 or v.NumSnowLayerNeg != 0:
        raise NotImplementedError(f"snowfall on the ground ({v.SnowfallGround} mm/s) is not ported")

    if v.FlagFrozenGround:
        dz1 = v.ThicknessSnowSoilLayer[0]
        v.SoilIce[0] = v.SoilIce[0] + (v.DewSoilSfcLiq - v.EvapSoilSfcLiq) * dt / (dz1 * 1000.0)
        v.DewSoilSfcLiq = 0.0
        v.EvapSoilSfcLiq = 0.0
        if v.SoilIce[0] < 0.0:
            v.SoilLiqWater[0] = v.SoilLiqWater[0] + v.SoilIce[0]
            v.SoilIce[0] = 0.0
        v.SoilMoisture[0] = v.SoilLiqWater[0] + v.SoilIce[0]
    # end-of-step total water as SoilWaterMain leaves it, up to the (not ported) change of liquid water
    v.SoilMoisture = [v.SoilLiqWater[k] + v.SoilIce[k] for k in range(n)]
    v.WaterToAtmosTotal = v.Transpiration + v.EvapCanopyNet + v.EvapGroundNet


def noahmp_main(v, p):
    """Ports noahmp/src/NoahmpMainMod.F90 (NoahmpMain) for our configuration: no irrigation (asserted: not cropland),
    no dynamic vegetation / crop biochemistry, balance checks skipped (diagnostic only)."""
    atmos.process_atmos_forcing(v, p)
    atmos.general_init(v, p)
    atmos.phenology_main(v, p)
    if v.FlagCropland:   # IrrigationPrepare can trigger sprinkler irrigation on cropland even with IRRIGATION_OPTION 0
        raise NotImplementedError("cropland (irrigation trigger) is not ported")
    canopy_water.canopy_water_intercept(v, p)
    canopy_water.precipitation_heat_advect(v, p)
    energy_main(v, p)
    water_main_subset(v, p)
