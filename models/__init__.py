import numpy as np

from model_builder.models.structure import mappings as mpg
from model_builder.models.structure import contacts as cts
from model_builder.models import potentials as ptl

'''
A Model consists of:
- a structure mapping
- a Hamiltonian
    - a set of parameterized interaction potentials
'''

class Model(object):
    """Model class """
    def __init__(self, topology, bead_repr="CA"):
        self.structure_mapping = mpg.assign_mapping(bead_repr, topology) 
        self.Hamiltonian = ptl.Hamiltonian()

    def describe(self):
        #TODO: What is the best description?
        pass 
    
    def map_traj(self, traj):
        return self.structure_mapping.map_traj(traj)

class StructureBasedModel(Model):

    def __init__(self, topology, bead_repr=None):
        Model.__init__(self, topology, bead_repr=bead_repr)
        self.Hamiltonian = ptl.StructureBasedHamiltonian()
         
    def set_reference(self, traj):
        self.ref_traj_aa = traj[0]
        self.ref_traj = self.structure_mapping.map_traj(traj[0])

    def add_sbm_potentials(self):
        self.Hamiltonian.add_sbm_potentials() 
        pass
