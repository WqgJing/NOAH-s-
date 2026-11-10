"""pyport: Python port of the part of Noah-MP (HRLDAS v5.2.1 + noah/patches/noahmp_sstar.patch) that our point
configuration uses. See README.md in this folder.

    from pyport.driver import PointModel
    m = PointModel("runs/pyref_VR")          # reads namelist.hrldas, NoahmpTable.TBL, setup, soil, forcing arrays
    times, out = m.run()                     # half-hourly values with LDASOUT names
"""
from .driver import PointModel, ldasout_record  # noqa: F401
