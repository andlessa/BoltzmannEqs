.. index:: Installation

Installation
============

The code is written in Python and requires Python version 2.6 or later (but not Python 3)
with the following *external* Python libraries:

 * `numpy <https://pypi.python.org/pypi/numpy>`_
 * `scipy <https://pypi.python.org/pypi/scipy>`_
 * `Assimulo <https://pypi.python.org/pypi/Assimulo>`_ >=2.9

In addition, some of the plotting functions may require:

 * `matplotlib <https://pypi.python.org/pypi/matplotlib>`_
 * `pyROOT <https://root.cern.ch/pyroot>`_

For a successful installation of Assilumo, `SUNDIALS <https://computation.llnl.gov/projects/sundials>`_ must also be installed.

A script for installing Sundials and Assimulo is provided: ::

   assimulo_installer.sh

The script downloads ths Sundials tarball and installs in a local folder (./sundials). It then installs Assimulo in the local user folder using pip. 
The installation can then be checked running: ::

   python check_installation.py

After the installation the user must add the path to the Sundials lib folder (./sundials/lib) to its enviroment variable (LD_LIBRARY_PATH).



