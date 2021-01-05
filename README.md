# Pore-Builder 
[![Build Status](https://dev.azure.com/rayamatsumoto/porebuilder/_apis/build/status/rmatsum836.Pore-Builder?branchName=refs%2Fpull%2F30%2Fmerge)](https://dev.azure.com/rayamatsumoto/porebuilder/_build/latest?definitionId=5&branchName=refs%2Fpull%2F30%2Fmerge)
[![codecov](https://codecov.io/gh/rmatsum836/Pore-Builder/branch/master/graph/badge.svg)](https://codecov.io/gh/rmatsum836/Pore-Builder) <br>

## Overview

`Pore-Builder` is an ![mBuild](http://mosdef-hub.github.io/mbuild/) recipe to build carbon slit pores and surfaces.
Recipes are designed to be add-ons, allowing users to develop routines to initialize specific chemistries using mBuild.

## Installation
To install `Pore-Builder`, its dependencies must first be installed.  It is recommended to do so through Conda by running the
following commands:
```
# Create conda environment
conda create -n pore37 python=3.7
# Add channels
conda config --add channels conda-forge
# Install packages
conda install --file requirements.txt
```
Once the dependencies have been installed, `Pore-Builder` can be installed from the package root directory: `pip install -e .`.
To check the package has been correctly installed, we can run Python and attempt to import the package via mBuild:
```
>>> import mbuild
>>> mbuild.recipes.GraphenePore()
/Users/raymatsumoto/anaconda3/envs/science37/lib/python3.7/site-packages/mbuild/lattice.py:626: UserWarning: Periodicity of non-rectangular lattices are not valid with default boxes. Only rectangular lattices are valid at this time.
  warn('Periodicity of non-rectangular lattices are not valid with '
<GraphenePore 2688 particles, periodicity: [3.9296     2.675      2.97774175], 0 bonds, id: 140234327632848>
>>> 
```

## Docker
A docker image is available on DockerHub to try out `Pore-Builder`.

## Publications using `Pore-Builder`
- https://www.tandfonline.com/doi/abs/10.1080/00268976.2020.1742938
- https://www.authorea.com/users/311036/articles/496195-open-source-molecular-modeling-software-in-chemical-engineering-focusing-on-the-molecular-simulation-design-framework?commit=bed95fe98da6c2ed7360a8ef4442e0b35ed82b63

![image](https://user-images.githubusercontent.com/25011342/68546370-838f9880-03a3-11ea-8db6-232c0d7a6dff.png)

