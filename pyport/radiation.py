"""Shortwave albedo / two-stream canopy radiation, absorbed and reflected shortwave, longwave emissivity.

RADIATIVE_TRANSFER_OPTION 3 (gap = 1 - FVEG), SNOW_ALBEDO_OPTION 1 (BATS), SNOW_COVER_OPTION 1 (Niu & Yang 2007).
"""
import math

import numpy as np

from .constants import TFRZ


def snow_cover_ground_niu07(v, p):
    """Ports noahmp/src/SnowCoverGroundNiu07Mod.F90 (SnowCoverGroundNiu07)."""
    v.SnowCoverFrac = 0.0
    if v.SnowDepth > 0.0:
        dens = v.SnowWaterEquiv / v.SnowDepth
        melt_fac = (dens / 100.0) ** p.SnowMeltFac
        v.SnowCoverFrac = math.tanh(v.SnowDepth / (p.SnowCoverFac * melt_fac))


def snow_aging_bats(v, p):
    """Ports noahmp/src/SnowAgingBatsMod.F90 (SnowAgingBats)."""
    if v.SnowWaterEquiv <= 0.0:
        v.SnowAgeNondim = 0.0
    else:
        time_fac = p.MainTimeStep / p.SnowAgeFacBats
        grow_exp = p.SnowGrowVapFacBats * (1.0 / TFRZ - 1.0 / v.TemperatureGrd)
        vap_eff = math.exp(grow_exp)
        frz_eff = math.exp(min(0.0, p.SnowGrowFrzFacBats * grow_exp))
        age_chg = time_fac * (vap_eff + frz_eff + p.SnowSootFacBats)
        fresh = max(0.0, v.SnowWaterEquiv - v.SnowWaterEquivPrev) / p.SnowMassFullCoverOld
        v.SnowAgeNondim = max(0.0, (v.SnowAgeNondim + age_chg) * (1.0 - fresh))
    v.SnowAgeFac = v.SnowAgeNondim / (v.SnowAgeNondim + 1.0)


def snow_albedo_bats(v, p):
    """Ports noahmp/src/SnowAlbedoBatsMod.F90 (SnowAlbedoBats)."""
    b = p.SolarZenithAdjBats
    zen = max((1.0 + 1.0 / b) / (1.0 + 2.0 * b * v.CosSolarZenithAngle) - 1.0 / b, 0.0)
    dif_vis = p.FreshSnoAlbVisBats * (1.0 - p.SnoAgeFacDifVisBats * v.SnowAgeFac)
    dif_nir = p.FreshSnoAlbNirBats * (1.0 - p.SnoAgeFacDifNirBats * v.SnowAgeFac)
    v.AlbedoSnowDif = np.array([dif_vis, dif_nir])
    v.AlbedoSnowDir = np.array([dif_vis + p.SzaFacDirVisBats * zen * (1.0 - dif_vis),
                                dif_nir + p.SzaFacDirNirBats * zen * (1.0 - dif_nir)])


def ground_albedo(v, p):
    """Ports noahmp/src/GroundAlbedoMod.F90 (GroundAlbedo), soil surface (SurfaceType 1)."""
    adj = max(0.11 - 0.40 * v.SoilMoisture[0], 0.0)
    soil_dir = np.minimum(p.AlbedoSoilSat + adj, p.AlbedoSoilDry)
    v.AlbedoSoilDir = soil_dir
    v.AlbedoSoilDif = soil_dir.copy()
    scf = v.SnowCoverFrac
    v.AlbedoGrdDir = soil_dir * (1.0 - scf) + v.AlbedoSnowDir * scf
    v.AlbedoGrdDif = v.AlbedoSoilDif * (1.0 - scf) + v.AlbedoSnowDif * scf


def canopy_radiation_two_stream(v, p, ib, idif):
    """Ports noahmp/src/CanopyRadiationTwoStreamMod.F90 (CanopyRadiationTwoStream) for band ib (0 vis, 1 nir) and
    idif (0 direct, 1 diffuse); gap fractions for RADIATIVE_TRANSFER_OPTION 3."""
    cosz = v.CosSolarZenithAngle
    vai = v.VegAreaIndEff
    if vai == 0.0:
        gap_dir, gap_dif = 1.0, 1.0
    else:
        if p.OptCanopyRadiationTransfer != 3:
            raise NotImplementedError("only RADIATIVE_TRANSFER_OPTION 3 is ported")
        gap_dir, gap_dif = 1.0 - v.VegFrac, 1.0 - v.VegFrac
    v.GapCanopyDir, v.GapCanopyDif = gap_dir, gap_dif

    cosz_t = max(0.001, cosz)
    xl = min(max(p.CanopyOrientIndex, -0.4), 0.6)
    if abs(xl) <= 0.01:
        xl = 0.01
    phi1 = 0.5 - 0.633 * xl - 0.330 * xl * xl
    phi2 = 0.877 * (1.0 - 2.0 * phi1)
    proj = phi1 + phi2 * cosz_t
    v.VegAreaProjDir = proj
    od_dir = proj / cosz_t
    od_dif = (1.0 - phi1 / phi2 * math.log((phi1 + phi2) / phi1)) / phi2
    rho, tau = v.ReflectanceVeg[ib], v.TransmittanceVeg[ib]
    omega_leaf = rho + tau
    tmp0 = proj + phi2 * cosz_t
    tmp1 = phi1 * cosz_t
    ssa = 0.5 * omega_leaf * proj / tmp0 * (1.0 - tmp1 / tmp0 * math.log((tmp1 + tmp0) / tmp1))
    up_dir_leaf = (1.0 + od_dif * od_dir) / (omega_leaf * od_dif * od_dir) * ssa
    up_dif_leaf = 0.5 * (rho + tau + (rho - tau) * ((1.0 + xl) / 2.0) ** 2) / omega_leaf

    # adjust omega, betad and betai for intercepted snow
    if v.TemperatureCanopy > TFRZ:
        t0, t1, t2 = omega_leaf, up_dir_leaf, up_dif_leaf
    else:
        fw = v.CanopyWetFrac
        oms = p.ScatterCoeffSnow[ib]
        t0 = (1.0 - fw) * omega_leaf + fw * oms
        t1 = ((1.0 - fw) * omega_leaf * up_dir_leaf + fw * oms * p.UpscatterCoeffSnowDir) / t0
        t2 = ((1.0 - fw) * omega_leaf * up_dif_leaf + fw * oms * p.UpscatterCoeffSnowDif) / t0
    om, betad, betai = t0, t1, t2

    # absorbed, reflected, transmitted fluxes per unit incoming radiation
    B = 1.0 - om + om * betai
    C = om * betai
    tmp0 = od_dif * od_dir
    D = tmp0 * om * betad
    F = tmp0 * om * (1.0 - betad)
    tmp1 = B * B - C * C
    H = math.sqrt(tmp1) / od_dif
    sigma = tmp0 * tmp0 - tmp1
    if abs(sigma) < 1.0e-6:
        sigma = math.copysign(1.0e-6, sigma)
    p1 = B + od_dif * H
    p2 = B - od_dif * H
    p3 = B + tmp0
    p4 = B - tmp0
    s1 = math.exp(-H * vai)
    s2 = math.exp(-od_dir * vai)
    alb_g = v.AlbedoGrdDir[ib] if idif == 0 else v.AlbedoGrdDif[ib]
    u1 = B - C / alb_g
    u2 = B - C * alb_g
    u3 = F + C * alb_g
    tmp2 = u1 - od_dif * H
    tmp3 = u1 + od_dif * H
    d1 = p1 * tmp2 / s1 - p2 * tmp3 * s1
    tmp4 = u2 + od_dif * H
    tmp5 = u2 - od_dif * H
    d2 = tmp4 / s1 - tmp5 * s1
    h1 = -D * p4 - C * F
    tmp6 = D - h1 * p3 / sigma
    tmp7 = (D - C - h1 / sigma * (u1 + tmp0)) * s2
    h2 = (tmp6 * tmp2 / s1 - p2 * tmp7) / d1
    h3 = -(tmp6 * tmp3 * s1 - p1 * tmp7) / d1
    h4 = -F * p3 - C * D
    tmp8 = h4 / sigma
    tmp9 = (u3 - tmp8 * (u2 - tmp0)) * s2
    h5 = -(tmp8 * tmp4 / s1 + tmp9) / d2
    h6 = (tmp8 * tmp5 * s1 + tmp9) / d2
    h7 = (C * tmp2) / (d1 * s1)
    h8 = (-C * tmp3 * s1) / d1
    h9 = tmp4 / (d2 * s1)
    h10 = (-tmp5 * s1) / d2

    # downward direct and diffuse fluxes below vegetation (Niu and Yang 2004)
    if idif == 0:
        trans_dir = s2 * (1.0 - gap_dir) + gap_dir
        trans_dif = (h4 * s2 / sigma + h5 * s1 + h6 / s1) * (1.0 - gap_dir)
        v.RadSwDirTranGrdDir[ib] = trans_dir
        v.RadSwDifTranGrdDir[ib] = trans_dif
        refl_tot = (h1 / sigma + h2 + h3) * (1.0 - gap_dir) + v.AlbedoGrdDir[ib] * gap_dir
        refl_can = (h1 / sigma + h2 + h3) * (1.0 - gap_dir)
        refl_grd = v.AlbedoGrdDir[ib] * gap_dir
        v.AlbedoSfcDir[ib] = refl_tot
        v.RadSwReflVegDir[ib] = refl_can
        v.RadSwReflGrdDir[ib] = refl_grd
        v.RadSwAbsVegDir[ib] = 1.0 - v.AlbedoSfcDir[ib] - (1.0 - v.AlbedoGrdDir[ib]) * v.RadSwDirTranGrdDir[ib] - \
            (1.0 - v.AlbedoGrdDif[ib]) * v.RadSwDifTranGrdDir[ib]
    else:
        trans_dir = 0.0
        trans_dif = (h9 * s1 + h10 / s1) * (1.0 - gap_dif) + gap_dif
        v.RadSwDirTranGrdDif[ib] = trans_dir
        v.RadSwDifTranGrdDif[ib] = trans_dif
        refl_tot = (h7 + h8) * (1.0 - gap_dif) + v.AlbedoGrdDif[ib] * gap_dif
        refl_can = (h7 + h8) * (1.0 - gap_dif) + v.AlbedoGrdDif[ib] * gap_dif
        v.AlbedoSfcDif[ib] = refl_tot
        v.RadSwReflVegDif[ib] = refl_can
        v.RadSwReflGrdDif[ib] = 0.0
        v.RadSwAbsVegDif[ib] = 1.0 - v.AlbedoSfcDif[ib] - (1.0 - v.AlbedoGrdDir[ib]) * v.RadSwDirTranGrdDif[ib] - \
            (1.0 - v.AlbedoGrdDif[ib]) * v.RadSwDifTranGrdDif[ib]


def surface_albedo(v, p):
    """Ports noahmp/src/SurfaceAlbedoMod.F90 (SurfaceAlbedo) with SNOW_ALBEDO_OPTION 1 (BATS)."""
    if p.OptSnowAlbedo != 1:
        raise NotImplementedError("only SNOW_ALBEDO_OPTION 1 (BATS) is ported")
    minthr = 1.0e-06
    v.GapBtwCanopy = 0.0
    v.GapInCanopy = 0.0
    v.VegAreaProjDir = 0.0
    v.ReflectanceVeg = np.zeros(2)
    v.TransmittanceVeg = np.zeros(2)
    v.CanopySunlitFrac = 0.0
    for name in ("AlbedoSfcDir", "AlbedoSfcDif", "AlbedoGrdDir", "AlbedoGrdDif", "AlbedoSnowDir", "AlbedoSnowDif",
                 "RadSwAbsVegDir", "RadSwAbsVegDif", "RadSwDirTranGrdDir", "RadSwDirTranGrdDif",
                 "RadSwDifTranGrdDir", "RadSwDifTranGrdDif", "RadSwReflVegDir", "RadSwReflVegDif",
                 "RadSwReflGrdDir", "RadSwReflGrdDif"):
        setattr(v, name, np.zeros(2))
    v.VegAreaIndEff = v.LeafAreaIndEff + v.StemAreaIndEff

    snow_aging_bats(v, p)

    cosz = v.CosSolarZenithAngle
    if cosz > 0:
        vai = v.VegAreaIndEff
        leaf_w = v.LeafAreaIndEff / max(vai, minthr)
        stem_w = v.StemAreaIndEff / max(vai, minthr)
        v.ReflectanceVeg = np.maximum(p.ReflectanceLeaf * leaf_w + p.ReflectanceStem * stem_w, minthr)
        v.TransmittanceVeg = np.maximum(p.TransmittanceLeaf * leaf_w + p.TransmittanceStem * stem_w, minthr)
        snow_albedo_bats(v, p)
        ground_albedo(v, p)
        for ib in range(p.NumSwRadBand):
            canopy_radiation_two_stream(v, p, ib, 0)
            canopy_radiation_two_stream(v, p, ib, 1)
        ext = v.VegAreaProjDir / cosz * math.sqrt(1.0 - v.ReflectanceVeg[0] - v.TransmittanceVeg[0])
        sun = (1.0 - math.exp(-ext * vai)) / max(ext * vai, minthr)
        v.CanopySunlitFrac = 0.0 if sun < 0.01 else sun

    v.CanopyShadeFrac = 1.0 - v.CanopySunlitFrac
    v.LeafAreaIndSunlit = v.LeafAreaIndEff * v.CanopySunlitFrac
    v.LeafAreaIndShade = v.LeafAreaIndEff * v.CanopyShadeFrac


def surface_radiation(v, p):
    """Ports noahmp/src/SurfaceRadiationMod.F90 (SurfaceRadiation), non-SNICAR path."""
    minthr = 1.0e-6
    abs_can_dir = np.zeros(2)
    abs_can_dif = np.zeros(2)
    v.RadSwAbsGrd = v.RadSwAbsVeg = v.RadSwAbsSfc = 0.0
    v.RadSwReflSfc = v.RadSwReflVeg = v.RadSwReflGrd = 0.0
    v.RadPhotoActAbsSunlit = v.RadPhotoActAbsShade = 0.0
    for ib in range(p.NumSwRadBand):
        abs_can_dir[ib] = v.RadSwDownDir[ib] * v.RadSwAbsVegDir[ib]
        abs_can_dif[ib] = v.RadSwDownDif[ib] * v.RadSwAbsVegDif[ib]
        v.RadSwAbsVeg += abs_can_dir[ib] + abs_can_dif[ib]
        v.RadSwAbsSfc += abs_can_dir[ib] + abs_can_dif[ib]
        tran_dir = v.RadSwDownDir[ib] * v.RadSwDirTranGrdDir[ib]
        tran_dif = v.RadSwDownDir[ib] * v.RadSwDifTranGrdDir[ib] + v.RadSwDownDif[ib] * v.RadSwDifTranGrdDif[ib]
        abs_grd = tran_dir * (1.0 - v.AlbedoGrdDir[ib]) + tran_dif * (1.0 - v.AlbedoGrdDif[ib])
        v.RadSwAbsGrd += abs_grd
        v.RadSwAbsSfc += abs_grd

    lai_frac = v.LeafAreaIndEff / max(v.VegAreaIndEff, minthr)
    if v.CanopySunlitFrac > 0.0:
        v.RadPhotoActAbsSunlit = (abs_can_dir[0] + v.CanopySunlitFrac * abs_can_dif[0]) * lai_frac / \
            max(v.LeafAreaIndSunlit, minthr)
        v.RadPhotoActAbsShade = (v.CanopyShadeFrac * abs_can_dif[0]) * lai_frac / max(v.LeafAreaIndShade, minthr)
    else:
        v.RadPhotoActAbsSunlit = 0.0
        v.RadPhotoActAbsShade = (abs_can_dir[0] + abs_can_dif[0]) * lai_frac / max(v.LeafAreaIndShade, minthr)

    d, f = v.RadSwDownDir, v.RadSwDownDif
    v.RadSwReflSfc = (v.AlbedoSfcDir[0] * d[0] + v.AlbedoSfcDif[0] * f[0]) + \
                     (v.AlbedoSfcDir[1] * d[1] + v.AlbedoSfcDif[1] * f[1])
    v.RadSwReflVeg = v.RadSwReflVegDir[0] * d[0] + v.RadSwReflVegDif[0] * f[0] + \
        v.RadSwReflVegDir[1] * d[1] + v.RadSwReflVegDif[1] * f[1]
    v.RadSwReflGrd = v.RadSwReflGrdDir[0] * d[0] + v.RadSwReflGrdDif[0] * f[0] + \
        v.RadSwReflGrdDir[1] * d[1] + v.RadSwReflGrdDif[1] * f[1]


def surface_emissivity(v, p):
    """Ports noahmp/src/SurfaceEmissivityMod.F90 (SurfaceEmissivity), soil point (IndicatorIceSfc = 0)."""
    ev = 1.0 - math.exp(-(v.LeafAreaIndEff + v.StemAreaIndEff) / 1.0)
    eg = p.EmissivitySoilLake[0] * (1.0 - v.SnowCoverFrac) + p.EmissivitySnow * v.SnowCoverFrac
    fveg = v.VegFrac
    v.EmissivityVeg, v.EmissivityGrd = ev, eg
    v.EmissivitySfc = fveg * (eg * (1 - ev) + ev + ev * (1 - ev) * (1 - eg)) + (1 - fveg) * eg
