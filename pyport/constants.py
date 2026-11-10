"""Physical constants of Noah-MP.

Ports noahmp/src/ConstantDefineMod.F90. The Fortran build is single precision (kind_noahmp = 4, utility/Machine.F90),
so every constant is stored here as the float32 value the Fortran code actually uses (e.g. 273.16 -> 273.1600036...),
then used in float64 arithmetic.
"""
import numpy as np


def f32(x):
    """Round a Python number to the nearest float32 value (what a Fortran `real` literal or variable holds)."""
    return float(np.float32(x))


GRAV = f32(9.80616)          # ConstGravityAcc [m/s2]
SB = f32(5.67e-08)           # ConstStefanBoltzmann [W/m2/K4]
VKC = f32(0.40)              # ConstVonKarman
TFRZ = f32(273.16)           # ConstFreezePoint [K]
HSUB = f32(2.8440e06)        # ConstLatHeatSublim [J/kg]
HVAP = f32(2.5104e06)        # ConstLatHeatEvap [J/kg]
HFUS = f32(0.3336e06)        # ConstLatHeatFusion [J/kg]
CWAT = f32(4.188e06)         # ConstHeatCapacWater [J/m3/K]
CICE = f32(2.094e06)         # ConstHeatCapacIce [J/m3/K]
CPAIR = f32(1004.64)         # ConstHeatCapacAir [J/kg/K]
TKWAT = f32(0.57)            # ConstThermConductWater [W/m/K]
TKICE = f32(2.2)             # ConstThermConductIce [W/m/K]
TKAIR = f32(0.023)           # ConstThermConductAir [W/m/K]
TKQTZ = f32(7.7)             # ConstThermConductQuartz [W/m/K]
TKO = f32(2.0)               # ConstThermConductSoilOth [W/m/K]
RAIR = f32(287.04)           # ConstGasDryAir [J/kg/K]
RW = f32(461.269)            # ConstGasWaterVapor [J/kg/K]
DENH2O = f32(1000.0)         # ConstDensityWater [kg/m3]
DENICE = f32(917.0)          # ConstDensityIce [kg/m3]
PI = f32(3.14159265)         # ConstPI

UNDEFINED = -9999.0          # undefined_real (utility/Machine.F90)
