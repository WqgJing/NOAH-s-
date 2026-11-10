"""Turbulent and radiative energy balance of the vegetated and bare tiles (SURFACE_DRAG_OPTION 1, Monin-Obukhov;
CANOPY_STOMATAL_RESISTANCE_OPTION 1, Ball-Berry; TEMP_TIME_SCHEME_OPTION 1 or 3)."""
import math

from .constants import f32, CPAIR, SB, VKC, GRAV, TFRZ, CWAT, CICE, DENH2O, DENICE

# saturation vapour pressure polynomial coefficients (VaporPressureSaturationMod.F90), float32 like the Fortran parameters
_A = [f32(x) for x in (6.107799961, 4.436518521e-01, 1.428945805e-02, 2.650648471e-04, 3.031240396e-06,
                       2.034080948e-08, 6.136820929e-11)]
_B = [f32(x) for x in (6.109177956, 5.034698970e-01, 1.886013408e-02, 4.176223716e-04, 5.824720280e-06,
                       4.838803174e-08, 1.838826904e-10)]
_C = [f32(x) for x in (4.438099984e-01, 2.857002636e-02, 7.938054040e-04, 1.215215065e-05, 1.036561403e-07,
                       3.532421810e-10, -7.090244804e-13)]
_D = [f32(x) for x in (5.030305237e-01, 3.773255020e-02, 1.267995369e-03, 2.477563108e-05, 3.005693132e-07,
                       2.158542548e-09, 7.131097725e-12)]


def _poly(c, t):
    return 100.0 * (c[0] + t * (c[1] + t * (c[2] + t * (c[3] + t * (c[4] + t * (c[5] + t * c[6]))))))


def vapor_pressure_saturation(t):
    """Ports noahmp/src/VaporPressureSaturationMod.F90 (VaporPressureSaturation). t in deg C.
    Returns (esat water, esat ice, desat/dT water, desat/dT ice) [Pa, Pa/K]."""
    return _poly(_A, t), _poly(_B, t), _poly(_C, t), _poly(_D, t)


def _temp_unit_conv(t):
    """Statement function TempUnitConv: Kelvin -> deg C limited to [-50, 50]."""
    return min(50.0, max(-50.0, t - TFRZ))


def _esat(t_kelvin):
    tc = _temp_unit_conv(t_kelvin)
    esw, esi, dsw, dsi = vapor_pressure_saturation(tc)
    return (esw, dsw) if tc > 0.0 else (esi, dsi)


def _psi_mh(zeta):
    """Stability functions used by the MOST routines (unstable: Paulson 1970; stable: -5 zeta)."""
    if zeta < 0.0:
        x = (1.0 - 16.0 * zeta) ** 0.25
        t2 = math.log((1.0 + x * x) / 2.0)
        t3 = math.log((1.0 + x) / 2.0)
        return 2.0 * t3 + t2 - 2.0 * math.atan(x) + 1.5707963, 2 * t2
    return -5.0 * zeta, -5.0 * zeta


# ---------------------------------------------------------------------------------------------------------------------
# vegetated tile
# ---------------------------------------------------------------------------------------------------------------------
def resistance_above_canopy_most(v, p, it, sh_tmp, sgn):
    """Ports noahmp/src/ResistanceAboveCanopyMostMod.F90 (ResistanceAboveCanopyMOST). `it` is the 1-based iteration
    index; returns the updated sign-change counter MoStabParaSgn."""
    mpe = 1.0e-6
    mozold = v.MoStabParaAbvCan
    zref, d0 = v.RefHeightAboveGrd, v.ZeroPlaneDispSfc
    if zref <= d0:
        raise RuntimeError("RefHeightAboveGrd <= ZeroPlaneDispSfc")
    tmpcm = math.log((zref - d0) / v.RoughLenMomSfc)
    tmpch = math.log((zref - d0) / v.RoughLenShCanopy)
    tmpcm2 = math.log((2.0 + v.RoughLenMomSfc) / v.RoughLenMomSfc)
    tmpch2 = math.log((2.0 + v.RoughLenShCanopy) / v.RoughLenShCanopy)
    if it == 1:
        v.FrictionVelVeg = 0.0
        v.MoStabParaAbvCan = 0.0
        v.MoLengthAbvCan = 0.0
        v.MoStabParaVeg2m = 0.0
    else:
        tvir = (1.0 + 0.61 * v.SpecHumidityRefHeight) * v.TemperatureAirRefHeight
        tmp1 = VKC * (GRAV / tvir) * sh_tmp / (v.DensityAirRefHeight * CPAIR)
        if abs(tmp1) <= mpe:
            tmp1 = mpe
        v.MoLengthAbvCan = -1.0 * v.FrictionVelVeg ** 3 / tmp1
        v.MoStabParaAbvCan = min((zref - d0) / v.MoLengthAbvCan, 1.0)
        v.MoStabParaVeg2m = min((2.0 + v.RoughLenShCanopy) / v.MoLengthAbvCan, 1.0)
    if mozold * v.MoStabParaAbvCan < 0.0:
        sgn += 1
    if sgn >= 2:
        v.MoStabParaAbvCan = 0.0
        v.MoStabCorrMomAbvCan = v.MoStabCorrShAbvCan = 0.0
        v.MoStabParaVeg2m = 0.0
        v.MoStabCorrMomVeg2m = v.MoStabCorrShVeg2m = 0.0
    fm, fh = _psi_mh(v.MoStabParaAbvCan)
    if v.MoStabParaAbvCan < 0.0:
        fm2, fh2 = _psi_mh(v.MoStabParaVeg2m)
    else:
        fm2 = fh2 = -5.0 * v.MoStabParaVeg2m
    if it == 1:
        v.MoStabCorrMomAbvCan, v.MoStabCorrShAbvCan = fm, fh
        v.MoStabCorrMomVeg2m, v.MoStabCorrShVeg2m = fm2, fh2
    else:
        v.MoStabCorrMomAbvCan = 0.5 * (v.MoStabCorrMomAbvCan + fm)
        v.MoStabCorrShAbvCan = 0.5 * (v.MoStabCorrShAbvCan + fh)
        v.MoStabCorrMomVeg2m = 0.5 * (v.MoStabCorrMomVeg2m + fm2)
        v.MoStabCorrShVeg2m = 0.5 * (v.MoStabCorrShVeg2m + fh2)
    v.MoStabCorrShAbvCan = min(v.MoStabCorrShAbvCan, 0.9 * tmpch)
    v.MoStabCorrMomAbvCan = min(v.MoStabCorrMomAbvCan, 0.9 * tmpcm)
    v.MoStabCorrShVeg2m = min(v.MoStabCorrShVeg2m, 0.9 * tmpch2)
    v.MoStabCorrMomVeg2m = min(v.MoStabCorrMomVeg2m, 0.9 * tmpcm2)
    cmfm = tmpcm - v.MoStabCorrMomAbvCan
    chfh = tmpch - v.MoStabCorrShAbvCan
    cm2fm2 = tmpcm2 - v.MoStabCorrMomVeg2m
    ch2fh2 = tmpch2 - v.MoStabCorrShVeg2m
    if abs(cmfm) <= mpe:
        cmfm = mpe
    if abs(chfh) <= mpe:
        chfh = mpe
    if abs(cm2fm2) <= mpe:
        cm2fm2 = mpe
    if abs(ch2fh2) <= mpe:
        ch2fh2 = mpe
    v.ExchCoeffMomAbvCan = VKC * VKC / (cmfm * cmfm)
    v.ExchCoeffShAbvCan = VKC * VKC / (cmfm * chfh)
    v.FrictionVelVeg = v.WindSpdRefHeight * math.sqrt(v.ExchCoeffMomAbvCan)
    v.ExchCoeffSh2mVegMo = VKC * v.FrictionVelVeg / ch2fh2
    v.ResistanceMomAbvCan = max(1.0, 1.0 / (v.ExchCoeffMomAbvCan * v.WindSpdRefHeight))
    v.ResistanceShAbvCan = max(1.0, 1.0 / (v.ExchCoeffShAbvCan * v.WindSpdRefHeight))
    v.ResistanceLhAbvCan = v.ResistanceShAbvCan
    return sgn


def resistance_leaf_to_ground(v, p, it, vai, sh_grd_tmp):
    """Ports noahmp/src/ResistanceLeafToGroundMod.F90 (ResistanceLeafToGround): under-canopy aerodynamic resistance
    and bulk leaf boundary-layer resistance."""
    mpe = 1.0e-6
    v.MoStabParaUndCan = 0.0
    v.MoLengthUndCan = 0.0
    if it > 1:
        tmp1 = VKC * (GRAV / v.TemperatureCanopyAir) * sh_grd_tmp / (v.DensityAirRefHeight * CPAIR)
        if abs(tmp1) <= mpe:
            tmp1 = mpe
        v.MoLengthUndCan = -1.0 * v.FrictionVelVeg ** 3 / tmp1
        v.MoStabParaUndCan = min((v.ZeroPlaneDispSfc - v.RoughLenMomGrd) / v.MoLengthUndCan, 1.0)
    if v.MoStabParaUndCan < 0.0:
        fhg = (1.0 - 15.0 * v.MoStabParaUndCan) ** (-0.25)
    else:
        fhg = 1.0 + 4.7 * v.MoStabParaUndCan
    if it == 1:
        v.MoStabCorrShUndCan = fhg
    else:
        v.MoStabCorrShUndCan = 0.5 * (v.MoStabCorrShUndCan + fhg)
    hc = v.CanopyHeight
    kext = (p.CanopyWindExtFac * vai * hc * v.MoStabCorrShUndCan) ** 0.5
    v.WindExtCoeffCanopy = kext
    tmp1 = math.exp(-kext * v.RoughLenShVegGrd / hc)
    tmp2 = math.exp(-kext * (v.RoughLenShCanopy + v.ZeroPlaneDispSfc) / hc)
    tmprah2 = hc * math.exp(kext) / kext * (tmp1 - tmp2)
    kh = max(VKC * v.FrictionVelVeg * (hc - v.ZeroPlaneDispSfc), mpe)
    v.ResistanceMomUndCan = 0.0
    v.ResistanceShUndCan = tmprah2 / kh
    v.ResistanceLhUndCan = v.ResistanceShUndCan
    tmprb = kext * 50.0 / (1.0 - math.exp(-kext / 2.0))
    rb = tmprb * math.sqrt(p.LeafDimLength / v.WindSpdCanopyTop)
    v.ResistanceLeafBoundary = min(max(rb, 5.0), 50.0)


def resistance_canopy_stomata_ball_berry(v, p, shade):
    """Ports noahmp/src/ResistanceCanopyStomataBallBerryMod.F90 (ResistanceCanopyStomataBallBerry); shade = 0 sunlit,
    1 shaded. BTRAN (SoilTranspFacAcc) multiplies Vcmax."""
    mpe = 1.0e-6
    f1 = lambda ab, bc: ab ** ((bc - 25.0) / 10.0)
    f2 = lambda ab: 1.0 + math.exp((-2.2e05 + 710.0 * (ab + 273.16)) / (8.314 * (ab + 273.16)))
    P, ta = v.PressureAirRefHeight, v.TemperatureAirRefHeight
    cf = P / (8.314 * ta) * 1.0e06
    rs = 1.0 / p.ConductanceLeafMin * cf
    psn = 0.0
    par = (v.RadPhotoActAbsSunlit if shade == 0 else v.RadPhotoActAbsShade) / max(v.VegFrac, 1.0e-6)
    if par > 0.0:
        c3 = p.PhotosynPathC3
        fnf = min(v.NitrogenConcFoliage / max(mpe, p.NitrogenConcFoliageMax), 1.0)
        tc = v.TemperatureCanopy - TFRZ
        ppf = 4.6 * par
        j = ppf * p.QuantumEfficiency25C
        kc = p.Co2MmConst25C * f1(p.Co2MmConstQ10, tc)
        ko = p.O2MmConst25C * f1(p.O2MmConstQ10, tc)
        awc = kc * (1.0 + v.PressureAtmosO2 / ko)
        cp = 0.5 * kc / ko * v.PressureAtmosO2 * 0.21
        vcmx = p.CarboxylRateMax25C / f2(tc) * fnf * v.SoilTranspFacAcc * f1(p.CarboxylRateMaxQ10, tc)
        co2 = v.PressureAtmosCO2
        ci = 0.7 * co2 * c3 + 0.4 * co2 * (1.0 - c3)
        rlb = v.ResistanceLeafBoundary / cf
        esat = v.VapPresSatCanopy
        cea = max(0.25 * esat * c3 + 0.40 * esat * (1.0 - c3), min(v.PressureVaporCanAir, esat))
        for _ in range(3):
            wj = max(ci - cp, 0.0) * j / (ci + 2.0 * cp) * c3 + j * (1.0 - c3)
            wc = max(ci - cp, 0.0) * vcmx / (ci + awc) * c3 + vcmx * (1.0 - c3)
            we = 0.5 * vcmx * c3 + 4000.0 * vcmx * ci / P * (1.0 - c3)
            psn = min(wj, wc, we) * v.IndexGrowSeason
            cs = max(co2 - 1.37 * rlb * P * psn, mpe)
            a = p.SlopeConductToPhotosyn * psn * P * cea / (cs * esat) + p.ConductanceLeafMin
            b = (p.SlopeConductToPhotosyn * psn * P / cs + p.ConductanceLeafMin) * rlb - 1.0
            c = -rlb
            if b >= 0.0:
                q = -0.5 * (b + math.sqrt(b * b - 4.0 * a * c))
            else:
                q = -0.5 * (b - math.sqrt(b * b - 4.0 * a * c))
            r1, r2 = q / a, c / q
            rs = max(r1, r2)
            ci = max(cs - psn * P * 1.65 * rs, 0.0)
        rs = rs * cf
    if shade == 0:
        v.ResistanceStomataSunlit, v.PhotosynLeafSunlit = rs, psn
    else:
        v.ResistanceStomataShade, v.PhotosynLeafShade = rs, psn


def surface_energy_flux_vegetated(v, p):
    """Ports noahmp/src/SurfaceEnergyFluxVegetatedMod.F90 (SurfaceEnergyFluxVegetated): canopy temperature by
    Newton-type iteration with M-O stability (up to 20 iterations), then the under-canopy ground temperature
    (5 iterations), 2-m diagnostics and output conductances."""
    if p.OptSurfaceDrag != 1:
        raise NotImplementedError("only SURFACE_DRAG_OPTION 1 is ported")
    if p.OptStomataResistance != 1:
        raise NotImplementedError("only CANOPY_STOMATAL_RESISTANCE_OPTION 1 is ported")
    dt = p.MainTimeStep
    rho, ta, P = v.DensityAirRefHeight, v.TemperatureAirRefHeight, v.PressureAirRefHeight
    fveg = v.VegFrac
    ev, eg = v.EmissivityVeg, v.EmissivityGrd
    lwdn = v.RadLwDownRefHeight

    last_iter = 0
    v.FrictionVelVeg = 0.1
    v.MoStabParaAbvCan = 0.0
    sgn = 0
    v.MoStabCorrShVeg2m = 0.0
    sh_grd_tmp = 0.0
    sh_can_tmp = 0.0
    vai = min(6.0, v.VegAreaIndEff)
    lai_sun = min(6.0, v.LeafAreaIndSunlit)
    lai_sha = min(6.0, v.LeafAreaIndShade)

    v.VapPresSatGrdVeg, _ = _esat(v.TemperatureGrdVeg)

    v.CanopyHeight = p.HeightCanopyTop
    hc = v.CanopyHeight
    v.WindSpdCanopyTop = v.WindSpdRefHeight * math.log((hc - v.ZeroPlaneDispSfc + v.RoughLenMomSfc) /
                                                       v.RoughLenMomSfc) / \
        math.log(v.RefHeightAboveGrd / v.RoughLenMomSfc)
    if (hc - v.ZeroPlaneDispSfc) <= 0.0:
        raise RuntimeError("CanopyHeight <= ZeroPlaneDispSfc")

    lw_air = -ev * (1.0 + (1.0 - ev) * (1.0 - eg)) * lwdn - ev * eg * SB * v.TemperatureGrdVeg ** 4
    lw_can = (2.0 - ev * (1.0 - eg)) * ev * SB

    for it in range(1, 21):
        v.RoughLenShCanopy = v.RoughLenMomSfc
        v.RoughLenShVegGrd = v.RoughLenMomGrd
        sgn = resistance_above_canopy_most(v, p, it, sh_can_tmp, sgn)
        resistance_leaf_to_ground(v, p, it, vai, sh_grd_tmp)

        v.VapPresSatCanopy, v.VapPresSatCanTempD = _esat(v.TemperatureCanopy)

        if it == 1:
            resistance_canopy_stomata_ball_berry(v, p, 0)
            resistance_canopy_stomata_ball_berry(v, p, 1)

        ch_abv = 1.0 / v.ResistanceShAbvCan
        ch_leaf = 2.0 * vai / v.ResistanceLeafBoundary
        grd_coeff = 1.0 / v.ResistanceShUndCan
        tot = ch_abv + ch_leaf + grd_coeff
        temp_shgh = (ta * ch_abv + v.TemperatureGrdVeg * grd_coeff) / tot
        sh_frac = ch_leaf / tot
        sh_coeff = (1.0 - sh_frac) * rho * CPAIR * ch_leaf

        v.ExchCoeffLhAbvCan = 1.0 / v.ResistanceLhAbvCan
        rb = v.ResistanceLeafBoundary
        pot_evap = fveg * rho * CPAIR * vai / rb * (v.VapPresSatCanopy - v.PressureVaporCanAir) / v.PsychConstCanopy
        if pot_evap > 0.0 and v.CanopyWetFrac > 0.0:
            store = v.CanopyLiqWater if v.TemperatureCanopy > TFRZ else v.CanopyIce
            v.ExchCoeffLhEvap = min(v.CanopyWetFrac, store * v.LatHeatVapCanopy / dt / pot_evap) * vai / rb
        else:
            v.ExchCoeffLhEvap = v.CanopyWetFrac * vai / rb
        v.ExchCoeffLhTransp = (1.0 - v.CanopyWetFrac) * (lai_sun / (rb + v.ResistanceStomataSunlit) +
                                                         lai_sha / (rb + v.ResistanceStomataShade))
        v.ExchCoeffLhUndCan = 1.0 / (v.ResistanceLhUndCan + v.ResistanceGrdEvap)
        tot = v.ExchCoeffLhAbvCan + v.ExchCoeffLhEvap + v.ExchCoeffLhTransp + v.ExchCoeffLhUndCan
        vp_lh = (v.PressureVaporRefHeight * v.ExchCoeffLhAbvCan + v.VapPresSatGrdVeg * v.ExchCoeffLhUndCan) / tot
        et_frac = (v.ExchCoeffLhEvap + v.ExchCoeffLhTransp) / tot
        lh_coeff = (1.0 - et_frac) * v.ExchCoeffLhEvap * rho * CPAIR / v.PsychConstCanopy
        tr_coeff = (1.0 - et_frac) * v.ExchCoeffLhTransp * rho * CPAIR / v.PsychConstCanopy

        tv = v.TemperatureCanopy
        v.TemperatureCanopyAir = temp_shgh + sh_frac * tv
        v.PressureVaporCanAir = vp_lh + et_frac * v.VapPresSatCanopy
        v.RadLwNetCanopy = fveg * (lw_air + lw_can * tv ** 4)
        v.HeatSensibleCanopy = fveg * rho * CPAIR * ch_leaf * (tv - v.TemperatureCanopyAir)
        v.HeatLatentCanEvap = fveg * rho * CPAIR * v.ExchCoeffLhEvap * \
            (v.VapPresSatCanopy - v.PressureVaporCanAir) / v.PsychConstCanopy
        v.HeatLatentCanTransp = fveg * rho * CPAIR * v.ExchCoeffLhTransp * \
            (v.VapPresSatCanopy - v.PressureVaporCanAir) / v.PsychConstCanopy
        if tv > TFRZ:
            v.HeatLatentCanEvap = min(v.CanopyLiqWater * v.LatHeatVapCanopy / dt, v.HeatLatentCanEvap)
        else:
            v.HeatLatentCanEvap = min(v.CanopyIce * v.LatHeatVapCanopy / dt, v.HeatLatentCanEvap)
        heat_cap = p.HeatCapacCanFac * vai * CWAT + v.CanopyLiqWater * CWAT / DENH2O + v.CanopyIce * CICE / DENICE
        resid = v.RadSwAbsVeg - v.RadLwNetCanopy - v.HeatSensibleCanopy - v.HeatLatentCanEvap - \
            v.HeatLatentCanTransp + v.HeatPrecipAdvCanopy
        tot_coeff = fveg * (4.0 * lw_can * tv ** 3 + sh_coeff + (lh_coeff + tr_coeff) * v.VapPresSatCanTempD +
                            heat_cap / dt)
        dtv = resid / tot_coeff
        v.RadLwNetCanopy += fveg * 4.0 * lw_can * tv ** 3 * dtv
        v.HeatSensibleCanopy += fveg * sh_coeff * dtv
        v.HeatLatentCanEvap += fveg * lh_coeff * v.VapPresSatCanTempD * dtv
        v.HeatLatentCanTransp += fveg * tr_coeff * v.VapPresSatCanTempD * dtv
        v.HeatCanStorageChg = fveg * heat_cap / dt * dtv
        v.TemperatureCanopy = tv + dtv

        sh_can_tmp = rho * CPAIR * (v.TemperatureCanopyAir - ta) / v.ResistanceShAbvCan
        sh_grd_tmp = rho * CPAIR * (v.TemperatureGrdVeg - v.TemperatureCanopyAir) / v.ResistanceShUndCan

        v.SpecHumiditySfc = (0.622 * v.PressureVaporCanAir) / (P - 0.378 * v.PressureVaporCanAir)
        if last_iter == 1:
            break
        if it >= 5 and abs(dtv) <= 0.01 and last_iter == 0:
            last_iter = 1
    v.NumIterVeg = it

    # under-canopy fluxes and ground temperature
    lw_air = -eg * (1.0 - ev) * lwdn - eg * ev * SB * v.TemperatureCanopy ** 4
    lw_can = eg * SB
    sh_coeff = rho * CPAIR / v.ResistanceShUndCan
    lh_coeff = rho * CPAIR / (v.PsychConstGrd * (v.ResistanceLhUndCan + v.ResistanceGrdEvap))
    grd_coeff = 2.0 * v.ThermConductSoilSnow[0] / v.ThicknessSnowSoilLayer[0]
    t1 = v.TemperatureSoilSnow[0]
    for _ in range(5):
        tg = v.TemperatureGrdVeg
        v.VapPresSatGrdVeg, v.VapPresSatGrdVegTempD = _esat(tg)
        v.RadLwNetVegGrd = lw_can * tg ** 4 + lw_air
        v.HeatSensibleVegGrd = sh_coeff * (tg - v.TemperatureCanopyAir)
        v.HeatLatentVegGrd = lh_coeff * (v.VapPresSatGrdVeg * v.RelHumidityGrd - v.PressureVaporCanAir)
        v.HeatGroundVegGrd = grd_coeff * (tg - t1)
        resid = v.RadSwAbsGrd - v.RadLwNetVegGrd - v.HeatSensibleVegGrd - v.HeatLatentVegGrd - \
            v.HeatGroundVegGrd + v.HeatPrecipAdvVegGrd
        tot_coeff = 4.0 * lw_can * tg ** 3 + sh_coeff + lh_coeff * v.VapPresSatGrdVegTempD + grd_coeff
        dtg = resid / tot_coeff
        v.RadLwNetVegGrd += 4.0 * lw_can * tg ** 3 * dtg
        v.HeatSensibleVegGrd += sh_coeff * dtg
        v.HeatLatentVegGrd += lh_coeff * v.VapPresSatGrdVegTempD * dtg
        v.HeatGroundVegGrd += grd_coeff * dtg
        v.TemperatureGrdVeg = tg + dtg

    if p.OptSnowSoilTempTime in (1, 3):
        if (v.SnowDepth > 0.05) and (v.TemperatureGrdVeg > TFRZ):
            raise NotImplementedError("snow on the ground (SnowDepth > 0.05 m) is not ported")

    v.WindStressEwVeg = -rho * v.ExchCoeffMomAbvCan * v.WindSpdRefHeight * v.WindEastwardRefHeight
    v.WindStressNsVeg = -rho * v.ExchCoeffMomAbvCan * v.WindSpdRefHeight * v.WindNorthwardRefHeight

    # 2-m temperature and humidity over vegetation
    v.ExchCoeffSh2mVeg = v.FrictionVelVeg * VKC / (math.log((2.0 + v.RoughLenShCanopy) / v.RoughLenShCanopy) -
                                                   v.MoStabCorrShVeg2m)
    if v.ExchCoeffSh2mVeg < 1.0e-5:
        v.TemperatureAir2mVeg = v.TemperatureCanopyAir
        v.SpecHumidity2mVeg = v.SpecHumiditySfc
    else:
        v.TemperatureAir2mVeg = v.TemperatureCanopyAir - (v.HeatSensibleVegGrd + v.HeatSensibleCanopy / fveg) / \
            (rho * CPAIR) * 1.0 / v.ExchCoeffSh2mVeg
        v.SpecHumidity2mVeg = v.SpecHumiditySfc - ((v.HeatLatentCanEvap + v.HeatLatentCanTransp) / fveg +
                                                   v.HeatLatentVegGrd) / \
            (v.LatHeatVapCanopy * rho) * 1.0 / v.ExchCoeffSh2mVeg

    v.ExchCoeffShAbvCan = ch_abv
    v.ExchCoeffShLeaf = ch_leaf
    v.ExchCoeffShUndCan = 1.0 / v.ResistanceShUndCan


# ---------------------------------------------------------------------------------------------------------------------
# bare tile
# ---------------------------------------------------------------------------------------------------------------------
def resistance_bare_ground_most(v, p, it, sh_tmp, sgn):
    """Ports noahmp/src/ResistanceBareGroundMostMod.F90 (ResistanceBareGroundMOST); returns MoStabParaSgn."""
    mpe = 1.0e-6
    mozold = v.MoStabParaBare
    zref, d0 = v.RefHeightAboveGrd, v.ZeroPlaneDispGrd
    if zref <= d0:
        raise RuntimeError("RefHeightAboveGrd <= ZeroPlaneDispGrd")
    tmpcm = math.log((zref - d0) / v.RoughLenMomGrd)
    tmpch = math.log((zref - d0) / v.RoughLenShBareGrd)
    tmpcm2 = math.log((2.0 + v.RoughLenMomGrd) / v.RoughLenMomGrd)
    tmpch2 = math.log((2.0 + v.RoughLenShBareGrd) / v.RoughLenShBareGrd)
    if it == 1:
        v.FrictionVelBare = 0.0
        v.MoStabParaBare = 0.0
        v.MoLengthBare = 0.0
        v.MoStabParaBare2m = 0.0
    else:
        tvir = (1.0 + 0.61 * v.SpecHumidityRefHeight) * v.TemperatureAirRefHeight
        tmp1 = VKC * (GRAV / tvir) * sh_tmp / (v.DensityAirRefHeight * CPAIR)
        if abs(tmp1) <= mpe:
            tmp1 = mpe
        v.MoLengthBare = -1.0 * v.FrictionVelBare ** 3 / tmp1
        v.MoStabParaBare = min((zref - d0) / v.MoLengthBare, 1.0)
        v.MoStabParaBare2m = min((2.0 + v.RoughLenShBareGrd) / v.MoLengthBare, 1.0)
    if mozold * v.MoStabParaBare < 0.0:
        sgn += 1
    if sgn >= 2:
        v.MoStabParaBare = 0.0
        v.MoStabCorrMomBare = v.MoStabCorrShBare = 0.0
        v.MoStabParaBare2m = 0.0
        v.MoStabCorrMomBare2m = v.MoStabCorrShBare2m = 0.0
    fm, fh = _psi_mh(v.MoStabParaBare)
    if v.MoStabParaBare < 0.0:
        fm2, fh2 = _psi_mh(v.MoStabParaBare2m)
    else:
        fm2 = fh2 = -5.0 * v.MoStabParaBare2m
    if it == 1:
        v.MoStabCorrMomBare, v.MoStabCorrShBare = fm, fh
        v.MoStabCorrMomBare2m, v.MoStabCorrShBare2m = fm2, fh2
    else:
        v.MoStabCorrMomBare = 0.5 * (v.MoStabCorrMomBare + fm)
        v.MoStabCorrShBare = 0.5 * (v.MoStabCorrShBare + fh)
        v.MoStabCorrMomBare2m = 0.5 * (v.MoStabCorrMomBare2m + fm2)
        v.MoStabCorrShBare2m = 0.5 * (v.MoStabCorrShBare2m + fh2)
    v.MoStabCorrShBare = min(v.MoStabCorrShBare, 0.9 * tmpch)
    v.MoStabCorrMomBare = min(v.MoStabCorrMomBare, 0.9 * tmpcm)
    v.MoStabCorrShBare2m = min(v.MoStabCorrShBare2m, 0.9 * tmpch2)
    v.MoStabCorrMomBare2m = min(v.MoStabCorrMomBare2m, 0.9 * tmpcm2)
    cmfm = tmpcm - v.MoStabCorrMomBare
    chfh = tmpch - v.MoStabCorrShBare
    cm2fm2 = tmpcm2 - v.MoStabCorrMomBare2m
    ch2fh2 = tmpch2 - v.MoStabCorrShBare2m
    if abs(cmfm) <= mpe:
        cmfm = mpe
    if abs(chfh) <= mpe:
        chfh = mpe
    if abs(cm2fm2) <= mpe:
        cm2fm2 = mpe
    if abs(ch2fh2) <= mpe:
        ch2fh2 = mpe
    v.ExchCoeffMomBare = VKC * VKC / (cmfm * cmfm)
    v.ExchCoeffShBare = VKC * VKC / (cmfm * chfh)
    v.FrictionVelBare = v.WindSpdRefHeight * math.sqrt(v.ExchCoeffMomBare)
    v.ExchCoeffSh2mBareMo = VKC * v.FrictionVelBare / ch2fh2
    v.ResistanceMomBareGrd = max(1.0, 1.0 / (v.ExchCoeffMomBare * v.WindSpdRefHeight))
    v.ResistanceShBareGrd = max(1.0, 1.0 / (v.ExchCoeffShBare * v.WindSpdRefHeight))
    v.ResistanceLhBareGrd = v.ResistanceShBareGrd
    return sgn


def surface_energy_flux_bare_ground(v, p):
    """Ports noahmp/src/SurfaceEnergyFluxBareGroundMod.F90 (SurfaceEnergyFluxBareGround): bare-ground temperature
    (5 iterations with M-O stability), 2-m diagnostics."""
    rho, ta, P = v.DensityAirRefHeight, v.TemperatureAirRefHeight, v.PressureAirRefHeight
    eg = v.EmissivityGrd
    lwdn = v.RadLwDownRefHeight
    v.MoStabParaBare = 0.0
    sgn = 0
    v.MoStabCorrShBare2m = 0.0
    v.FrictionVelBare = 0.1
    sh_tmp = 0.0
    lw_coeff = eg * SB
    grd_coeff = 2.0 * v.ThermConductSoilSnow[0] / v.ThicknessSnowSoilLayer[0]
    t1 = v.TemperatureSoilSnow[0]
    for it in range(1, 6):
        v.RoughLenShBareGrd = v.RoughLenMomGrd
        sgn = resistance_bare_ground_most(v, p, it, sh_tmp, sgn)
        ch_mom_tmp = 1.0 / v.ResistanceMomBareGrd
        ch_sh_tmp = 1.0 / v.ResistanceShBareGrd
        tg = v.TemperatureGrdBare
        v.VapPresSatGrdBare, v.VapPresSatGrdBareTempD = _esat(tg)
        sh_coeff = rho * CPAIR / v.ResistanceShBareGrd
        lh_coeff = rho * CPAIR / v.PsychConstGrd / (v.ResistanceGrdEvap + v.ResistanceLhBareGrd)
        v.RadLwNetBareGrd = lw_coeff * tg ** 4 - eg * lwdn
        v.HeatSensibleBareGrd = sh_coeff * (tg - ta)
        v.HeatLatentBareGrd = lh_coeff * (v.VapPresSatGrdBare * v.RelHumidityGrd - v.PressureVaporRefHeight)
        v.HeatGroundBareGrd = grd_coeff * (tg - t1)
        resid = v.RadSwAbsGrd - v.RadLwNetBareGrd - v.HeatSensibleBareGrd - v.HeatLatentBareGrd - \
            v.HeatGroundBareGrd + v.HeatPrecipAdvBareGrd
        tot_coeff = 4.0 * lw_coeff * tg ** 3 + sh_coeff + lh_coeff * v.VapPresSatGrdBareTempD + grd_coeff
        dtg = resid / tot_coeff
        v.RadLwNetBareGrd += 4.0 * lw_coeff * tg ** 3 * dtg
        v.HeatSensibleBareGrd += sh_coeff * dtg
        v.HeatLatentBareGrd += lh_coeff * v.VapPresSatGrdBareTempD * dtg
        v.HeatGroundBareGrd += grd_coeff * dtg
        v.TemperatureGrdBare = tg + dtg
        sh_tmp = sh_coeff * (v.TemperatureGrdBare - ta)
        v.VapPresSatGrdBare, _ = _esat(v.TemperatureGrdBare)
        e_sfc = v.VapPresSatGrdBare * v.RelHumidityGrd
        v.SpecHumiditySfc = 0.622 * e_sfc / (P - 0.378 * e_sfc)
        v.MoistureFluxSfc = (v.SpecHumiditySfc - v.SpecHumidityRefHeight) * lh_coeff * v.PsychConstGrd / CPAIR

    if p.OptSnowSoilTempTime in (1, 3):
        if (v.SnowDepth > 0.05) and (v.TemperatureGrdBare > TFRZ):
            raise NotImplementedError("snow on the ground (SnowDepth > 0.05 m) is not ported")

    v.WindStressEwBare = -rho * v.ExchCoeffMomBare * v.WindSpdRefHeight * v.WindEastwardRefHeight
    v.WindStressNsBare = -rho * v.ExchCoeffMomBare * v.WindSpdRefHeight * v.WindNorthwardRefHeight

    v.ExchCoeffSh2mBare = v.FrictionVelBare * VKC / (math.log((2.0 + v.RoughLenShBareGrd) / v.RoughLenShBareGrd) -
                                                     v.MoStabCorrShBare2m)
    if v.ExchCoeffSh2mBare < 1.0e-5:
        v.TemperatureAir2mBare = v.TemperatureGrdBare
        v.SpecHumidity2mBare = v.SpecHumiditySfc
    else:
        v.TemperatureAir2mBare = v.TemperatureGrdBare - v.HeatSensibleBareGrd / (rho * CPAIR) * 1.0 / \
            v.ExchCoeffSh2mBare
        v.SpecHumidity2mBare = v.SpecHumiditySfc - v.HeatLatentBareGrd / (v.LatHeatVapGrd * rho) * \
            (1.0 / v.ExchCoeffSh2mBare + v.ResistanceGrdEvap)
    v.ExchCoeffShBare = ch_sh_tmp
    v.ExchCoeffMomBareDiag = ch_mom_tmp
