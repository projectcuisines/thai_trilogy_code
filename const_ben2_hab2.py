# -*- coding: utf-8 -*-
"""Physical constants for the Ben2 and Hab2/Hab2* experiments of the THAI project."""
# ATMOSPHERIC COMPOSITION
ptot = 1e5  # sufrace pressure [Pa]

# Ben2 and Hab2/Hab2* specific constants
co2 = 1.0  # mixing ratio of CO2 [bar]
n2 = 1 - co2  # mixing ratio of N2 [bar]
c_p = 846  # specific heat at constant volume [J K-1]
rcp = 0.189  # ratio molecular gas contant by specific gas constant
rgas = 188.92  # specific gas constant [J kg-1 K]

# General atmospheric constants
rvapor = 461.52  # specific gas constant of water vapor [J kg-1 K]
mw_dryair = 28.0134 * n2 + 44.01 * co2  # atmospheric molecular weight [g/mol]
mw_vapor = 18.01528  # water molecular weight [g/mol]
mw_ratio = mw_dryair / mw_vapor  # ratio of mw_dryair by mw_vapor
latent_heat = 2_501_000  # latent heat of vaporization [J kg-1]
mgas_constant = 8.31446261815324  # molecular gas constant [J kg-1 mol-1]

# PLANETS PROPERTIES
planet_dfactor = 0.920  # TRAPPIST-1e diameter [Earth diameter unit]
planet_gfactor = 0.930  # TRAPPIST-1e gravity [Earth diameter unit]

earth_d = 6371 * 2  # Earth diameter [km]
earth_g = 9.81  # Earth gravity [m s-2]
gplanet = planet_gfactor * earth_g  # TRAPPIST-1e gravity [m s-2]
dplanet = planet_dfactor * earth_d  # TRAPPIST-1e diameter [km]
rplanet = dplanet / 2  # TRAPPIST-1e radius [km] (km needed for PSG)
rplanet_m = rplanet * 1e3

# ORBITAL INPUTS
semiaxis = 28.17e-3  # TRAPPIST-1e semi major axis [AU]
period = 6.099615  # TRAPPIST-1e orbital period [day]
transit_duration = 345  # TRAPPIST-1e transit duration [s]

number_transits = 1  # Number of transits integrated by PSG
transit_sum = transit_duration * number_transits  # total integration time [s]
noisetime = 0.2  # exposure time [s]
noiseframes = transit_sum / noisetime  # number of exposure

# STAR INPUTS
orig_temp = 2566  # TRAPPIST-1 temperature [K] from Agol et al 2020
sun_r = 695700  # Sun radius [km]
radius = 0.1192  # TRAPPIST-1 radius [Sun radius unit] from Kahne et al., 2018
