# SVJ flat analysis

Python-style package to analyze ntuples created with TreeMaker.

To be replaced with `coffea.NanoEvents` at some point (see https://github.com/CoffeaTeam/coffea/blob/d5126ca0f5ebc00522d4fc453cb7635fb73b3d9c/coffea/nanoevents/schemas.py#L290 ), but some technical debt has been incurred.


## Installation

Recommended setup is a conda environment with at least the following packages:

```
conda install -c conda-forge root
conda install -c conda-forge notebook
pip install uproot
python -m pip install matplotlib
conda install -c conda-forge tqdm 
pip install mplhep
pip install coffea
```

Then install this package:

```
git clone git@github.com:boostedsvj/svjflatanalysis.git
pip install -e svjflatanalysis
```
