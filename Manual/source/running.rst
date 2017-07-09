.. _running:

Running the code
================

A basic example of how to run the code is provided
by Example.py. If all dependencies have been
successfully installed, this code can be run as: ::

   Example.py -p parameters.ini -o output.dat -P

The resulting output will be written to the output file output.dat.
If the flag -P is set, a simple plot showing the evolution of the energy
densities will also be displayed.
The main arguments are described below:

**usage**: nonThermalDM.py [-h] [-p PARAMETERFILE] [-o OUTPUTFILE] [-P]


*optional arguments*:
  -h, --help            show this help message and exit
  -p PARAMETERFILE, --parameterFile PARAMETERFILE
                        name of parameter file, where most options are defined
  -o OUTPUTFILE, --outputFile OUTPUTFILE
                        name of output file (optional argument). If not
                        define, no output will be saved
  -P, --plotResult      show simple plot for the evolution of densities



Model Definitions
-----------------

For each new BSM species, the user must supply the main functions required
by the Boltzmann equations. These functions are supplied when defining the
new state in the main file (using the `Component class <../_build/html/code.html#component.Component>`_ ) and
correspond to:

  * **decays** : specify the BRs and total width. Can be temperature-dependent.
  * **mass** : specify the mass. Can be temperature-dependent.
  * **sigv** : specify the thermally averaged annihilation cross-section. Can be temperature-dependent.

For the example provided these functions are given in *modelDefinitions.py*.

  
Furthermore, for each species, the user must specify if its spin set by the *dof* property
:math:`dof = \pm 2*spin + 1` (where the plus sign is for bosons and minus for fermions).


.. _parameterFile:


The Parameters File
-------------------

The basic models parameters are set in the parameters file.
These are model dependent and should be adapted for each individual input model.
For the simple models assumed in nonThermalDM.py, the following parameters are defined:

* *model parameters*:

  * **yCoupling** (float): Mediator-DM-SM coupling
  * **mMediator** (float): Mediator mass (in GeV)
    perfoming :ref:`mass compression <massComp>`. *Only used if doCompress = True*
  * **mDM** (float): DM mass (in GeV)
  * **TRH** (float): reheat temperature (in GeV)
  * **TF** (float): Final temperature for evolution (should be after mediator decay)


The Output
----------

If an output file is defined, all the output will be written to the file.
This output includes the input parameters used, a summary of the main cosmological
observables (relic densities, effective number of additional neutrinos,...).
Furthermore, a table containing the values for the number and energy densities
of each species defined as a function of the temperature (or scale factor) is also printed.

If no output file is defined, only the summary information will be written to the file.
An example of a typical output is the following:

.. literalinclude:: output_example.dat
   :lines: 1-24

...

