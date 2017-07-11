# Boltzmann Code - PQMSSM #

This branch focus on an implementation of the supersymmetric axion model (PQMSSM).
It is significantly more complex than the version in master, since
it requires as input the whole MSSM (through an SLHA file), as well as the axion, saxion
and axino properties.
This version still is under construction and needs to be updated.

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


## Running BoltzmannCode ##

A basic example of how to run the code is provided
by Example.py. If all dependencies have been
successfully installed, this code can be run as: ::

   Example.py -p parameters.ini -o output.dat -P

The resulting output will be written to the output file output.dat.
If the flag -P is set, a simple plot showing the evolution of the energy
densities will also be displayed.

## Manual ##

For a detailed documentation (currently under construction), [read the docs](http://boltzmanncode.readthedocs.io)

### Contact ###

For questions, comments and bug reports, contact Andre Lessa at andre.lessa@ufabc.edu.br
