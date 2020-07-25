[![Python 3.8](https://img.shields.io/badge/python-3.8-blue.svg?logo=python&logoColor=white)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-black.svg)](LICENSE)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

:warning: **Work in progress** :warning:

Python code for analysing THAI data.

Contributions are welcome.

## Setting up
1. Install Miniconda or Anaconda.
2. In the command line, navigate to this folder and type
```bash
conda env create --file environment.yml
```

## Description of the output parameters

| Variable | Units | UM name |
|----------|-------|------------------|
| OLR clear| W m-2 | toa_outgoing_longwave_flux_assuming_clear_sky |
| OLR cloudy| W m-2 | toa_outgoing_longwave_flux |
| ASR clear| W m-2 | toa_incoming_shortwave_flux - toa_outgoing_shortwave_flux_assuming_clear_sky |
| ASR cloudy| W m-2 | toa_incoming_shortwave_flux - toa_outgoing_shortwave_flux |
| Surface temperature map | K | surface_temperature |
| Downward total SW flux at the surface | W m-2 | surface_downwelling_shortwave_flux_in_air |
| Net LW flux at the surface | W m-2 | surface_net_downward_longwave_flux |
| Open ocean fraction| km2 | OPEN OCEAN FRACTION or STASH_m01s03i395 |
| Total cloud liquid water path | kg m-2 | atmosphere_cloud_liquid_water_content |
| Total cloud ice water path | kg m-2 | atmosphere_cloud_ice_content |
| Total vapor column | kg m-2 | TOTAL COLUMN Q (WATER VAPOUR PATH) or STASH_m01s30i461 |
| Atmospheric temperature | K | air_temperature |
| U wind speed | m s-1 | eastward_wind |
| V wind speed | m s-1 | northward_wind |
| W wind speed | m s-1 | upward_air_velocity |
| 3D Heating rates SW | K s-1 | tendency_of_air_temperature_due_to_shortwave_heating |
| 3D Heating rates LW | K s-1 | tendency_of_air_temperature_due_to_longwave_heating |
| 3D Specific humidity | kg kg-1 | specific_humidity |
| 3D relative humidity | % | relative_humidity |
| Total CF | 0-1 | BULK CLOUD FRACTION IN EACH LAYER or STASH_m01s00i266 |
| Liquid CF | 0-1 | LIQUID CLOUD FRACTION IN EACH LAYER or STASH_m01s00i267 |
| Ice CF | 0-1 | FROZEN CLOUD FRACTION IN EACH LAYER or STASH_m01s00i268 |
| Total cloud MMR | kg kg-1 | mass_fraction_of_cloud_liquid_water_in_air + mass_fraction_of_cloud_ice_in_air |
| Liquid cloud MMR | kg kg-1 | mass_fraction_of_cloud_liquid_water_in_air |
| Ice cloud MMR | kg kg-1 | mass_fraction_of_cloud_ice_in_air |
| Pressure at rho (u/v) levels | Pa | air_pressure |
| Pressure at theta levels | Pa | air_pressure |
| OLR in each band | W m-2 | UPWARD LW FLUX ON LEVELS AND BANDS |
