import numpy as np

from model_builder.models.structure import mappings as mpg
from model_builder.models.structure import contacts as cts
from model_builder.models import potentials as ptl

'''
A Model consists of:
- a structure mapping
- a set of potentials
    - a set of interaction parameters
'''

class Model(object):

    def __init__(self, topology=None, traj=None):
        if (topology is None) and (traj is not None):
            topology = traj.top
        self.structure_mapping = mpg.CalphaMapping(topology)  
        self.potentials = ptl.Hamiltonian()

    def set_reference(self, traj):
        self.ref_traj_aa = traj[0]
        self.ref_traj = self.structure_mapping.map_traj(traj[0])

    def describe(self):
        pass 

    def add_sbm_contacts(self):
        self.potentials.add_sbm_contacts(self)

    def Vij(self):
        return [ interaction.Vij for interaction in self.pairV ]

    def map_traj(self, traj):
        self.structure_mapping.map_traj(traj)

