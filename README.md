# Boltzmann Code - FIMP #

This branch focus on a simple non-thermal (FIMP) DM implementation.

## Installation ##

The code is written in Python and requires Python version 2.6 or later (but not Python 3)
with the following *external* Python libraries:

 * [numpy](https://pypi.python.org/pypi/numpy)
 * [scipy](https://pypi.python.org/pypi/scipy)
 * [Assimulo](https://pypi.python.org/pypi/Assimulo) >=2.9

In addition, some of the plotting functions may require:

 * [matplotlib](https://pypi.python.org/pypi/matplotlib)
 * [pyROOT](https://root.cern.ch/pyroot)

For a successful installation of Assilumo, [SUNDIALS](https://computation.llnl.gov/projects/sundials) must also be installed.

A script for installing Sundials and Assimulo is provided: 

```
assimulo_installer.sh
```

The script downloads the Sundials tarball and installs in a local folder (./sundials).
It then installs Assimulo in the local user folder using pip.
The installation can then be checked running:

```
python check_installation.py
```

After the installation the user must add the path to the Sundials lib folder (./sundials/lib)
to its enviroment variable (LD_LIBRARY_PATH).



## Running ##

A basic example of how to run the code is provided
by Example.py. If all dependencies have been
successfully installed, this code can be run as:


```
   nonThermalDM.py -p parameters.ini -o output.dat -P
```

The resulting output will be written to the output file output.dat.
If the flag -P is set, a simple plot showing the evolution of the energy
densities will also be displayed.

## Manual ##

For a detailed documentation (currently under construction), [read the docs](http://boltzmanncode.readthedocs.io)

### Contact ###

For questions, comments and bug reports, contact Andre Lessa at andre.lessa@ufabc.edu.br
