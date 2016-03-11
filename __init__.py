""" Package for generating Gromacs input files for coarse-grain Hamiltonians.

Description:
    
    This package generates the necessary Gromacs topology files for simulations
of coarse-grain models. The available models, which can be found in the models
submodule, are mostly Go-like models. This package is meant as a stand-alone.


Submodules:

models -- Holds classes of different coarse-grain models.


Example:

>>> ls
>>> 

"""

import models
import inputs
load_model = inputs.load_model
#import info_string
import make_model
