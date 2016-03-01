


'''
To Do:
The structure submodule should also handle adding cofactors to the system (e.g.
heme).
'''

import mappings
import viz_bonds
import contacts

from model_builder.models.structure.mappings import *

MAPPINGS = {"CA":CalphaMapping, "CACB":CalphaCbetaMapping}

def assign_mapping(code, topology):
    return MAPPINGS[code](topology)
