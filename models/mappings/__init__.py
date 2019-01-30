


'''
To Do:
The structure submodule should also handle adding cofactors to the system (e.g.
heme).
'''

from __future__ import absolute_import
from model_builder.models.mappings.calpha import CalphaMapping
from model_builder.models.mappings.calphacbeta import CalphaCbetaMapping
from model_builder.models.mappings.heavyatom import HeavyAtomMapping
from model_builder.models.mappings.awsem import AwsemMapping, AwsemBackboneMapping

def assign_mapping(code, topology):
    MAPPINGS = {"CA":CalphaMapping, "CACB":CalphaCbetaMapping,
                "All-Atom":HeavyAtomMapping, "AWSEM":AwsemMapping, 
                "AWSEM_backbone":AwsemBackboneMapping}
    return MAPPINGS[code](topology)
