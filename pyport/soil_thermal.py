"""Soil temperature (implicit thermal diffusion + tridiagonal solver) and soil water phase change, without snow layers.
TBOT_OPTION 2 (TMN at ZBOT), TEMP_TIME_SCHEME_OPTION 1 or 3, SUPERCOOLED_WATER_OPTION 1 (Niu & Yang 2006)."""
from .constants import HFUS, TFRZ, GRAV


def soil_snow_thermal_diffusion(v, p):
    """Ports noahmp/src/SoilSnowThermalDiffusionMod.F90 (SoilSnowThermalDiffusion). Returns the tridiagonal
    coefficients (A1, A2, A3, R) per unit time."""
    n = p.NumSoilLayer
    z = v.DepthSnowSoilLayer
    T, tk, hc = v.TemperatureSoilSnow, v.ThermConductSoilSnow, v.HeatCapacSoilSnow
    rad_pen = [0.0] * n                                   # RadSwPenetrateGrd (non-SNICAR: 0)
    hcpa, dinv, tgd, ee = [0.0] * n, [0.0] * n, [0.0] * n, [0.0] * n
    for k in range(n):
        if k == 0:
            hcpa[k] = -z[k] * hc[k]
            dz = -z[k + 1]
            dinv[k] = 2.0 / dz
            tgd[k] = 2.0 * (T[k] - T[k + 1]) / dz
            ee[k] = tk[k] * tgd[k] - v.HeatGroundTotMean - rad_pen[k]
        elif k < n - 1:
            hcpa[k] = (z[k - 1] - z[k]) * hc[k]
            dz = z[k - 1] - z[k + 1]
            dinv[k] = 2.0 / dz
            tgd[k] = 2.0 * (T[k] - T[k + 1]) / dz
            ee[k] = (tk[k] * tgd[k] - tk[k - 1] * tgd[k - 1]) - rad_pen[k]
        else:
            hcpa[k] = (z[k - 1] - z[k]) * hc[k]
            if p.OptSoilTemperatureBottom == 1:
                v.HeatFromSoilBot = 0.0
            elif p.OptSoilTemperatureBottom == 2:
                tgd[k] = (T[k] - v.TemperatureSoilBottom) / (0.5 * (z[k - 1] + z[k]) - v.DepthSoilTempBotToSno)
                v.HeatFromSoilBot = -tk[k] * tgd[k]
            ee[k] = (-v.HeatFromSoilBot - tk[k - 1] * tgd[k - 1]) - rad_pen[k]
    a1, a2, a3, r = [0.0] * n, [0.0] * n, [0.0] * n, [0.0] * n
    for k in range(n):
        if k == 0:
            a1[k] = 0.0
            a3[k] = -tk[k] * dinv[k] / hcpa[k]
            if p.OptSnowSoilTempTime in (1, 3):
                a2[k] = -a3[k]
            else:
                raise NotImplementedError("TEMP_TIME_SCHEME_OPTION 2 is not ported")
        elif k < n - 1:
            a1[k] = -tk[k - 1] * dinv[k - 1] / hcpa[k]
            a3[k] = -tk[k] * dinv[k] / hcpa[k]
            a2[k] = -(a1[k] + a3[k])
        else:
            a1[k] = -tk[k - 1] * dinv[k - 1] / hcpa[k]
            a3[k] = 0.0
            a2[k] = -(a1[k] + a3[k])
        r[k] = ee[k] / (-hcpa[k])
    return a1, a2, a3, r


def matrix_solver_tridiagonal(A, B, C, D):
    """Ports noahmp/src/MatrixSolverTriDiagonalMod.F90 (MatrixSolverTriDiagonal), top layer index 1. Returns the
    solution P (the Fortran P array)."""
    n = len(B)
    C = list(C)
    C[n - 1] = 0.0
    P, delta = [0.0] * n, [0.0] * n
    P[0] = -C[0] / B[0]
    delta[0] = D[0] / B[0]
    for k in range(1, n):
        fac = 1.0 / (B[k] + A[k] * P[k - 1])
        P[k] = -C[k] * fac
        delta[k] = (D[k] - A[k] * delta[k - 1]) * fac
    P[n - 1] = delta[n - 1]
    for kk in range(n - 2, -1, -1):
        P[kk] = P[kk] * P[kk + 1] + delta[kk]
    return P


def soil_snow_temperature_solver(v, p, dt, a1, a2, a3, r):
    """Ports noahmp/src/SoilSnowTemperatureSolverMod.F90 (SoilSnowTemperatureSolver)."""
    n = p.NumSoilLayer
    r = [x * dt for x in r]
    a1 = [x * dt for x in a1]
    a2 = [1.0 + x * dt for x in a2]
    a3 = [x * dt for x in a3]
    dT = matrix_solver_tridiagonal(a1, a2, a3, r)
    v.TemperatureSoilSnow = [v.TemperatureSoilSnow[k] + dT[k] for k in range(n)]


def soil_snow_temperature_main(v, p):
    """Ports noahmp/src/SoilSnowTemperatureMainMod.F90 (SoilSnowTemperatureMain), non-SNICAR (no solar penetration)."""
    v.DepthSoilTempBotToSno = p.DepthSoilTempBottom - v.SnowDepth
    a1, a2, a3, r = soil_snow_thermal_diffusion(v, p)
    soil_snow_temperature_solver(v, p, p.SoilTimeStep, a1, a2, a3, r)
    v.HeatFromSoilBot = v.HeatFromSoilBot * p.SoilTimeStep


def soil_water_supercool_niu06(p, k, t):
    """Ports noahmp/src/SoilWaterSupercoolNiu06Mod.F90 (SoilWaterSupercoolNiu06): max liquid water [m3/m3] at t < TFRZ."""
    pot = HFUS * (TFRZ - t) / (GRAV * t)
    return p.SoilMoistureSat[k] * (pot / p.SoilMatPotentialSat[k]) ** (-1.0 / p.SoilExpCoeffB[k])


def soil_snow_water_phase_change(v, p):
    """Ports noahmp/src/SoilSnowWaterPhaseChangeMod.F90 (SoilSnowWaterPhaseChange) for soil layers only (no snow
    layer, no snow without layer: the caller asserts SWE = 0)."""
    if p.OptSoilSupercoolWater != 1:
        raise NotImplementedError("only SUPERCOOLED_WATER_OPTION 1 is ported")
    if v.NumSnowLayerNeg != 0 or v.SnowWaterEquiv > 0.0:
        raise NotImplementedError("snow is not ported")
    n, dt = p.NumSoilLayer, p.MainTimeStep
    dz = v.ThicknessSnowSoilLayer
    T = list(v.TemperatureSoilSnow)
    liq = [v.SoilLiqWater[k] * dz[k] * 1000.0 for k in range(n)]
    ice = [(v.SoilMoisture[k] - v.SoilLiqWater[k]) * dz[k] * 1000.0 for k in range(n)]
    ice0, liq0 = list(ice), list(liq)
    tot0 = [ice[k] + liq[k] for k in range(n)]
    supercool = [0.0] * n
    for k in range(n):
        if T[k] < TFRZ:
            supercool[k] = soil_water_supercool_niu06(p, k, T[k]) * dz[k] * 1000.0
    idx = [0] * n
    for k in range(n):
        if ice[k] > 0.0 and T[k] >= TFRZ:
            idx[k] = 1
        if liq[k] > supercool[k] and T[k] < TFRZ:
            idx[k] = 2
    eres, wpc = [0.0] * n, [0.0] * n
    for k in range(n):
        if idx[k] > 0:
            eres[k] = (T[k] - TFRZ) / v.PhaseChgFacSoilSnow[k]
            T[k] = TFRZ
        if idx[k] == 1 and eres[k] < 0.0:
            eres[k], idx[k] = 0.0, 0
        if idx[k] == 2 and eres[k] > 0.0:
            eres[k], idx[k] = 0.0, 0
        wpc[k] = eres[k] * dt / HFUS
    v.HeatLhTotPhsChg = 0.0
    for k in range(n):
        if idx[k] > 0 and abs(eres[k]) > 0.0:
            left = 0.0
            if wpc[k] > 0.0:
                ice[k] = max(0.0, ice0[k] - wpc[k])
                left = eres[k] - HFUS * (ice0[k] - ice[k]) / dt
            elif wpc[k] < 0.0:
                if tot0[k] < supercool[k]:
                    ice[k] = 0.0
                else:
                    ice[k] = min(tot0[k] - supercool[k], ice0[k] - wpc[k])
                    ice[k] = max(ice[k], 0.0)
                left = eres[k] - HFUS * (ice0[k] - ice[k]) / dt
            liq[k] = max(0.0, tot0[k] - ice[k])
            if abs(left) > 0.0:
                T[k] = T[k] + v.PhaseChgFacSoilSnow[k] * left
            v.HeatLhTotPhsChg += HFUS * (ice0[k] - ice[k]) / dt
    v.TemperatureSoilSnow = T
    v.SoilLiqWater = [liq[k] / (1000.0 * dz[k]) for k in range(n)]
    v.SoilMoisture = [(liq[k] + ice[k]) / (1000.0 * dz[k]) for k in range(n)]
    v.IndexPhaseChange = idx
