[![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg?logo=python&logoColor=white)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-black.svg)](LICENSE)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Python code for analysing [THAI data](https://thai.emac.gsfc.nasa.gov/organization/thai).

Contributions are welcome.

## What's in this repo?
| File | Purpose |
|:-----|:--------|
| `names/` | Model-specific names of variables and coordinates |
| `THAI-Cloud-Distribution.ipynb` | Cloud maps and meridional profiles |
| `THAI-Global-Diag.ipynb` | Global diagnostics (tables)|
| `THAI-Meridional-Profiles.ipynb` | Meridional profiles of temperature, WVP, CWP and cloud fraction |
| `THAI-Radiation-Fluxes.ipynb` | Maps of radiation fluxes at TOA and the surface |
| `THAI-Rot-Div-Components.ipynb` | Helmholtz decomposition |
| `THAI-Save-Time-Mean.ipynb` | Average data in time and save in an intermediate dataset |
| `THAI-Save-Time-Series.ipynb` | Average selected variables spatially and save in an intermediate dataset |
| `THAI-Sfc-Temp-Hist.ipynb` | Histograms of surface temperature |
| `THAI-Substellar-Profiles.ipynb` | Profiles of temperature, lapse rate, radiative heating at the substellar point |
| `THAI-TL-Streamfunction.ipynb` | Mass streamfunction in tidally locked coordinates |
| `THAI-Time-Series.ipynb` | Time series |
| `THAI-Time-Spectrum.ipynb` | Time variability as a spectrum |
| `THAI-Zonal-Mean-Zonal-Wind.ipynb` | Zonally averaged zonal wind |
| `calc.py` | Calculation of various climate diagnostics |
| `commons.py` | Common project-specific objects |
| `const_ben1_hab1.py` | Constants for the Ben 1 and Hab 1 cases |
| `const_ben2_hab2.py` | Constants for the Ben 2 and Hab 2 cases |
| `grid.py` | Operations on GCMs' grids |
| `load_thai.py` | Functions to load and merge THAI datasets |
| `model_exocam.py` | Functions relevant only to ExoCAM |
| `model_lmdg.py` | Functions relevant only to LMD-G |
| `model_rocke3d.py` | Functions relevant only to ROCKE-3D |
| `model_um.py` | Functions relevant only to the UM |
| `mypaths.py` | Paths to files and directories relevant to this project |
| `paper.mplstyle` | Matplotlib stylesheet |
| `plot_func.py` | Plotting functions |

## Want to re-use the code? Here are the instructions.
### Set up Python
Skip the first two steps if you have Jupyter Lab with `nb_conda_kernels` installed already.
1. Install [Miniforge](https://github.com/conda-forge/miniforge).
2. Install necessary packages to the `base` environment. Make sure you are installing them from the `conda-forge` channel.
```bash
conda install jupyterlab nb_conda_kernels
```
3. Check out or download this repository.
4. In the command line, navigate to the downloaded folder and create a separate conda environment.
```bash
conda env create --file environment.yml
```

### Download THAI data
All GCM THAI data are permanently available for download [here](https://ckan.emac.gsfc.nasa.gov/organization/thai), with variables described for each dataset.
If you use those data in your own research, please cite the relevant paper and add the following statement:

> THAI data have been obtained from https://ckan.emac.gsfc.nasa.gov/organization/thai, a data repository of the Sellers Exoplanet Environments Collaboration (SEEC), which is funded in part by the NASA Planetary Science Divisions Internal Scientist Funding Model.

### Open the code
1. Start the Jupyter Lab, for example from the command line:
```bash
jupyter lab
```
2. Open noteboks in the `thai` environment and start coding.

### Pre-process data
If you want to re-draw figures from the THAI trilogy papers, we recommend that you create an intermediate dataset first by averaging the THAI data in time or in space.
To get a dataset of time-averaged variables, use the [THAI-Save-Time-Mean.ipynb](THAI-Save-Time-Mean.ipynb) notebook.
To get a dataset of selected time series, use the [THAI-Save-Time-Series.ipynb](THAI-Save-Time-Series.ipynb) notebook.
**Before running these notebooks, make sure the paths to THAI data are set up correctly in [mypaths.py](mypaths.py)!**
Also note that these two notebooks use [dask](https://dask.org) which uses multiple cores to speed up calculations. You can adapt the number of cores to be more suitable for your machine.
