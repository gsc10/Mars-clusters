# Meteoroid Fragmentation in the Martian Atmosphere and the Formation of Crater Clusters 

## Overview

This repository contains source code and data to accompany a manuscript submitted to JGR-Planets:

Meteoroid Fragmentation in the Martian Atmosphere and the Formation of Crater Clusters
by G. S. Collins, E. L. Newland, D. Schwarz, M. Coleman, S. McMullan, I. J. Daubar, Katarina Miljkovi\'c, Tanja Neidhart, Eleanor Sansom

The contents of the repository is as follows:

* `PaperFigures.ipynb` A Jupyter (Python 3.7) notebook documenting how all the figures in the paper were produced.
* `mars_clusters.py` The Python implementation of the Monte Carlo model used to generate all the synthetic crater clusters described and presented in the paper. Note that this script imports and calls the fragment-cloud package, which must be downloaded and installed separately (see below).
* `gtools` A Python module with supplementary functions used in both `mars_clusters.py` and `PaperFigures.ipynb`.
* `obs-data` The observational data of new impact craters and clusters on Mars from Neidhart et al. (submitted) and Daubar et al. (submitted) used in this work to compare with model results.
* `model-data` Synthetic crater cluster output data used to generate the figures in the paper. The model input parameters of specific models are documented in `record.csv` the summary statistics of each model is documented in `statistics.csv`. Output from individual Monte Carlo simulations are included as `output-<model_number>.csv`. `frag_summary-198.csv` is a file containing information about all surviving fragments for model 198.
* `crater-data` Model output for model 198 including individual locations of all crater in every cluster. 

## Installation

To run the Monte Carlo simulation you will need to download and install the fragment-cloud package. Brief installation instructions for that package are provided here for ease.

### Prerequisites

You need to have the following packages installed:

* `python3` v. 3.7 or later with `numpy` and the python development headers
* `C++17` compiler, e.g. `g++` or `clang`
* `cmake` v. 3.12 or later
* `boost` v. 1.70 or later, including the `python3`, `numpy`, and `unit_test_framework` components.

### Installation Instructions

1. Open command line and navigate to the `fragment-cloud/` folder.
2. Inside the `fragment-cloud/` folder, create a new folder called `debug/` or `release/` and `cd` into it.
3. Run the following command: `cmake ..` This will check all dependencies and create a file called `setup.py` inside the `fragment-cloud/` folder.
4. Navigate to the `fragment-cloud/` folder.

For the final step, please install the `fcm` software like a pip package: `python3 setup.py install`. Then the `fcm` package will be available in your default python path.

## Example usage

`mars_clusters.py` takes a number of command-line options to specify model input parameters and options.

Output files will be written to the current directory. A record of input parameters is appended to `record.csv`. The outcome of each impact scenario is stored in `output-<id>.csv`. A summary of the statistics of the Monte Carlo run is appended to `statistics.csv`.

To run a test saved with model ID "test" for just 10 scenarios:

```
>python3 mars_clusters --id test -i 10
```

For more informtion on input parameters and model options:

```
>python3 mars_clusters.py --help

Program for simulating distribution of small clusters/craters on Mars

optional arguments:
  -h, --help            show this help message and exit
  --id MODEL_ID, -I MODEL_ID
  --restart, -R
  --verbose, -V
  --craters, -C
  --scaling CSCALING, -S CSCALING
  --lift LIFT_COEF, -L LIFT_COEF
  --fragments, -F
  --velocity VELOCITY, -v VELOCITY
  --angle ANGLE, -t ANGLE
  --atmos ATMOS, -A ATMOS
  --break BREAKM, -B BREAKM
  --alpha ALPHA, -a ALPHA
  --minmass MINMASS, -m MINMASS
  --dlo DLO, -l DLO
  --dup DUP, -u DUP
  --med_strength MED_STRENGTH, -s MED_STRENGTH
  --width WIDTH, -w WIDTH
  --minrad CRAT_MIN, -c CRAT_MIN
  --ablation CAB, -e CAB
  --drag CD, -d CD
  --ss_disp SS_DISP, -x SS_DISP
  --fm_disp FM_DISP, -y FM_DISP
  --fvcmin FVC_MIN, -f FVC_MIN
  --fvcmax FVC_MAX, -g FVC_MAX
  --semin SE_MIN, -j SE_MIN
  --semax SE_MAX, -k SE_MAX
  --impacts IMPACTS, -i IMPACTS
  --numc NUMC, -n NUMC
```

