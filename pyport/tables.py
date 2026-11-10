"""Readers for the Noah-MP parameter table, the HRLDAS namelist, the setup file and the spatial soil file.

Ports what the Fortran reads at start-up for our configuration:
  * noahmp/drivers/hrldas/NoahmpReadTableMod.F90  (NoahmpTable.TBL: MODIS vegetation block, rad, global, general and
    STAS soil blocks; the USGS block is not read when MMINLU = MODIFIED_IGBP_MODIS_NOAH)
  * noahmp/drivers/hrldas/NoahmpReadNamelistMod.F90 (namelist.hrldas, &NOAHLSM_OFFLINE)
  * hrldas/IO_code/module_hrldas_netcdf_io.F: READLAND_HRLDAS, READINIT_HRLDAS, READVEG_HRLDAS (hrldas_setup.nc) and
    READ_3D_SOIL (soil_spatial.nc, SOIL_DATA_OPTION = 4)
All real values are rounded to float32, as the Fortran stores them (kind_noahmp = 4).
"""
import re

import numpy as np

from .constants import f32


# ---------------------------------------------------------------------------------------------------------------------
# Fortran namelist text parsing (only what NoahmpTable.TBL and namelist.hrldas use)
# ---------------------------------------------------------------------------------------------------------------------
def _strip_comment(line):
    """Remove a trailing Fortran '!' comment that is not inside quotes."""
    out, quote = [], None
    for ch in line:
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out)


def _parse_value(tok):
    tok = tok.strip()
    if tok[:1] in "'\"":
        return tok.strip("'\"")
    low = tok.lower()
    if low in (".true.", "t", ".t."):
        return True
    if low in (".false.", "f", ".f."):
        return False
    try:
        return int(tok)
    except ValueError:
        return float(low.replace("d", "e"))


def _split_values(text):
    """Split a namelist value list on commas / blanks, keeping quoted strings together."""
    vals = re.findall(r"'[^']*'|\"[^\"]*\"|[^,\s]+", text)
    return [_parse_value(v) for v in vals]


def read_namelist_block(path, block):
    """Return {NAME: [values]} for namelist group `&block` (case-insensitive names, upper-cased keys).

    Array-element assignments such as soil_thick_input(2) = 0.3 are stored under 'SOIL_THICK_INPUT(2)'.
    """
    text = open(path).read()
    m = re.search(rf"^\s*&{block}\b", text, re.M | re.I)
    if m is None:
        raise KeyError(f"namelist block &{block} not found in {path}")
    body = []
    for line in text[m.end():].splitlines():
        if re.match(r"^\s*/\s*$", line):
            break
        body.append(_strip_comment(line))
    body = " ".join(body)
    out = {}
    names = list(re.finditer(r"([A-Za-z_][A-Za-z0-9_]*(?:\(\s*\d+\s*\))?)\s*=", body))
    for i, nm in enumerate(names):
        end = names[i + 1].start() if i + 1 < len(names) else len(body)
        key = re.sub(r"\s+", "", nm.group(1)).upper()
        out[key] = _split_values(body[nm.end():end])
    return out


# ---------------------------------------------------------------------------------------------------------------------
# NoahmpTable.TBL
# ---------------------------------------------------------------------------------------------------------------------
def read_table(path, veg_type, soil_type_table, soil_color=4):
    """Read the parameters our configuration uses for one vegetation type (MODIS) and the STAS soil class.

    veg_type: IVGTYP (MODIS/IGBP index, 1-based); soil_type_table: ISLTYP (1-based, only used for the initial
    SMOIS cap and frozen-soil initialization in NoahmpInitMain and for QUARTZ before the spatial override);
    soil_color: fixed to 4 in ConfigVarInTransferMod.F90.
    Returns a dict of scalars / small arrays (float32-rounded).
    """
    veg = read_namelist_block(path, "noahmp_modis_parameters")
    rad = read_namelist_block(path, "noahmp_rad_parameters")
    glo = read_namelist_block(path, "noahmp_global_parameters")
    gen = read_namelist_block(path, "noahmp_general_parameters")
    cat = read_namelist_block(path, "noahmp_stas_soil_categories")
    if cat["SLTYPE"][0].strip() != "STAS":
        raise ValueError("only the STAS soil table is supported")
    soil = read_namelist_block(path, "noahmp_soil_stas_parameters")
    iv, isl, ic = veg_type - 1, soil_type_table - 1, soil_color - 1

    def v(name):
        return f32(veg[name][iv])

    T = {}
    # indices of special land types (MODIS)
    for name in ("ISURBAN", "ISWATER", "ISBARREN", "ISICE", "ISCROP", "EBLFOREST", "NATURAL", "URBTYPE_BEG"):
        T[name] = int(veg[name][0])
    for name in ("CH2OP", "DLEAF", "Z0MVT", "HVT", "HVB", "DEN", "RC", "MFSNO", "SCFFAC", "CBIOM", "XL", "CWPVT",
                 "C3PSN", "KC25", "AKC", "KO25", "AKO", "AVCMX", "BP", "MP", "QE25", "VCMX25", "FOLNMX", "TMIN", "RGL",
                 "RS", "HS", "TOPT", "RSMAX", "SLA"):
        T[name] = v(name)
    T["NROOT"] = int(veg["NROOT"][iv])
    T["RHOL"] = np.array([v("RHOL_VIS"), v("RHOL_NIR")])
    T["RHOS"] = np.array([v("RHOS_VIS"), v("RHOS_NIR")])
    T["TAUL"] = np.array([v("TAUL_VIS"), v("TAUL_NIR")])
    T["TAUS"] = np.array([v("TAUS_VIS"), v("TAUS_NIR")])
    # radiation block
    T["ALBSAT"] = np.array([f32(rad["ALBSAT_VIS"][ic]), f32(rad["ALBSAT_NIR"][ic])])
    T["ALBDRY"] = np.array([f32(rad["ALBDRY_VIS"][ic]), f32(rad["ALBDRY_NIR"][ic])])
    T["ALBLAK"] = np.array([f32(x) for x in rad["ALBLAK"][:2]])
    T["OMEGAS"] = np.array([f32(x) for x in rad["OMEGAS"][:2]])
    T["BETADS"] = f32(rad["BETADS"][0])
    T["BETAIS"] = f32(rad["BETAIS"][0])
    T["EG"] = np.array([f32(x) for x in rad["EG"][:2]])
    T["EICE"] = f32(rad["EICE"][0])
    # global block
    for name in ("CO2", "O2", "Z0SNO", "SNOW_EMIS", "SWEMX", "TAU0", "GRAIN_GROWTH", "EXTRA_GROWTH", "DIRT_SOOT",
                 "BATS_COSZ", "BATS_VIS_NEW", "BATS_NIR_NEW", "BATS_VIS_AGE", "BATS_NIR_AGE", "BATS_VIS_DIR",
                 "BATS_NIR_DIR", "RSURF_SNOW", "RSURF_EXP", "SNOWDEN_MAX", "PSIWLT", "Z0SOIL", "Z0LAKE"):
        T[name] = f32(glo[name][0])
    # general block
    T["CSOIL"] = f32(gen["CSOIL_DATA"][0])
    T["ZBOT"] = f32(gen["ZBOT_DATA"][0])
    T["CZIL"] = f32(gen["CZIL_DATA"][0])
    # STAS soil table for the setup soil class (initialization only; SOIL_DATA_OPTION 4 overrides in the physics)
    T["BEXP_TABLE"] = f32(soil["BB"][isl])
    T["SMCMAX_TABLE"] = f32(soil["MAXSMC"][isl])
    T["PSISAT_TABLE"] = f32(soil["SATPSI"][isl])
    T["QUARTZ_TABLE"] = f32(soil["QTZ"][isl])
    return T


# ---------------------------------------------------------------------------------------------------------------------
# namelist.hrldas
# ---------------------------------------------------------------------------------------------------------------------
def read_hrldas_namelist(path):
    """Read &NOAHLSM_OFFLINE with the defaults of NoahmpReadNamelistMod.F90 for keys that are absent."""
    raw = read_namelist_block(path, "NOAHLSM_OFFLINE")
    defaults = {
        "DYNAMIC_VEG_OPTION": 4, "CANOPY_STOMATAL_RESISTANCE_OPTION": 1, "BTR_OPTION": 1, "SURFACE_RUNOFF_OPTION": 3,
        "SUBSURFACE_RUNOFF_OPTION": 3, "SURFACE_DRAG_OPTION": 1, "SUPERCOOLED_WATER_OPTION": 1,
        "FROZEN_SOIL_OPTION": 1, "RADIATIVE_TRANSFER_OPTION": 3, "SNOW_ALBEDO_OPTION": 1,
        "SNOW_THERMAL_CONDUCTIVITY": 1, "SNOW_COMPACTION_OPTION": 2, "SNOW_COVER_OPTION": 1,
        "PCP_PARTITION_OPTION": 1, "TBOT_OPTION": 2, "TEMP_TIME_SCHEME_OPTION": 1, "GLACIER_OPTION": 1,
        "SURFACE_RESISTANCE_OPTION": 1, "SOIL_DATA_OPTION": 1, "PEDOTRANSFER_OPTION": 1, "CROP_OPTION": 0,
        "IRRIGATION_OPTION": 0, "IRRIGATION_METHOD": 0, "DVIC_INFILTRATION_OPTION": 1, "TILE_DRAINAGE_OPTION": 0,
        "WETLAND_OPTION": 0, "SOIL_TIMESTEP": 0.0, "ZLVL": 10.0, "SF_URBAN_PHYSICS": 0,
        "PRESCRIBE_SOIL_MOISTURE": 0, "LAIO_SSTAR": -9999.0, "LAIO_SWILT": -9999.0,
        "NUIS_GS_SCALE": 1.0, "NUIS_RSURF_SCALE": 1.0, "NUIS_Z0_SCALE": 1.0, "SPINUP_LOOPS": 0,
        "SKIP_FIRST_OUTPUT": False, "START_MIN": 0,
    }
    nl = dict(defaults)
    for k, val in raw.items():
        nl[k] = val[0] if len(val) == 1 else val
    nsoil = int(nl["NSOIL"])
    nl["SOIL_THICK_INPUT"] = np.array([f32(raw[f"SOIL_THICK_INPUT({k})"][0]) for k in range(1, nsoil + 1)])
    for k in ("ZLVL", "LAIO_SSTAR", "LAIO_SWILT", "NUIS_GS_SCALE", "NUIS_RSURF_SCALE", "NUIS_Z0_SCALE",
              "SOIL_TIMESTEP"):
        nl[k] = f32(nl[k])
    return nl


# ---------------------------------------------------------------------------------------------------------------------
# netCDF inputs
# ---------------------------------------------------------------------------------------------------------------------
def read_setup(path):
    """hrldas_setup.nc, single point: READLAND_HRLDAS + READINIT_HRLDAS + READVEG_HRLDAS (module_hrldas_netcdf_io.F).

    The soil levels of the setup file equal the namelist layers here (0.1/0.3/0.6/1.0 m), so init_interp is the
    identity (checked by the caller).
    """
    import netCDF4
    with netCDF4.Dataset(path) as nc:
        g = lambda n: np.asarray(nc[n][:], dtype=np.float64).ravel()
        S = {n: float(g(n)[0]) for n in ("XLAT", "XLONG", "TMN", "SHDMAX", "SHDMIN", "LAI", "SNOW", "CANWAT", "TSK")}
        S["SNODEP"] = float(g("SNODEP")[0]) if "SNODEP" in nc.variables else None
        S["IVGTYP"] = int(g("IVGTYP")[0])
        S["ISLTYP"] = int(g("ISLTYP")[0])
        S["XLAND"] = int(g("XLAND")[0])
        S["DZS"] = g("DZS")
        S["TSLB"] = g("TSLB")
        S["SMOIS"] = g("SMOIS")
        S["MMINLU"] = nc.getncattr("MMINLU")
    return S


def read_soil_spatial(path):
    """soil_spatial.nc (SOIL_DATA_OPTION = 4): READ_3D_SOIL. Returns per-layer float64 arrays (float32 values)."""
    import netCDF4
    with netCDF4.Dataset(path) as nc:
        return {n: np.asarray(nc[n][:], dtype=np.float64)[:, 0, 0]
                for n in ("bexp", "smcdry", "smcwlt", "smcref", "smcmax", "dksat", "dwsat", "psisat", "quartz")}
