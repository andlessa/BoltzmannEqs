#!/usr/bin/env python 


"""
Check the Assimulo installation
"""

import sys,os

assimulo_path = None

try:
    import assimulo
    assimulo_path = os.path.dirname(os.path.abspath(assimulo.__file__))
except:
    print "Assimulo installation failed."
    sys.exit()

try:
    from assimulo.solvers import sundials
except:
    if os.path.isfile(os.path.join(assimulo_path,"solvers/sundials.so")):
        print "It seems Assimulo can no longer find the SUNDIALS libraries. Check if it has been added to the path."
    else:
        print "It seems Assimulo could not find the SUNDIALS libraries during installation. Check if it SUNDIALS has been properly installed and added to the path."
    sys.exit()

print "INSTALLATION SUCCESSFUL"
