[![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg?logo=python&logoColor=white)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-black.svg)](LICENSE)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

:warning: **Work in progress** :warning:

Python code for analysing THAI data.

Contributions are welcome.

## Setting up Python
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

## Download THAI data
All GCM THAI data are permanently available for download [here](https://ckan.emac.gsfc.nasa.gov/organization/thai), with variables described for each dataset.
If you use those data in your own research, please cite the relevant paper and add the following statement:

> THAI data have been obtained from https://ckan.emac.gsfc.nasa.gov/organization/thai, a data repository of the Sellers Exoplanet Environments Collaboration (SEEC), which is funded in part by the NASA Planetary Science Divisions Internal Scientist Funding Model.

## Open the code
1. Start the Jupyter Lab, for example from the command line:
```bash
jupyter lab
```
2. Open noteboks in the `thai` environment and start coding.

## Process data
If you want to re-draw figures from the THAI trilogy, we recommend that you pre-process THAI data first by averaging the data in time or in space.
To get a dataset of time-averaged variables, execute the [THAI-Save-Time-Mean.ipynb](THAI-Save-Time-Mean.ipynb) notebook.
To get a dataset of selected time series, use the [THAI-Save-Time-Series.ipynb](THAI-Save-Time-Series.ipynb) notebook.
**Before running these notebooks, make sure the paths to THAI data are set up correctly in `mypaths.py`!**
