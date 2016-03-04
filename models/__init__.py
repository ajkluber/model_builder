'''
Classes
-------
Model

StructureBasedModel

'''

import numpy as np

import output

import structure.mappings as mpg
import structure.contacts as cts
import potentials as ptl


class Model(object):
    """Model class """
    def __init__(self, topology, bead_repr="CA"):
        self.mapping = mpg.assign_mapping(bead_repr, topology) 
        self.Hamiltonian = ptl.Hamiltonian()

    def describe(self):
        #TODO: What is the best description?
        pass 
    
    def map_traj(self, traj):
        return self.mapping.map_traj(traj)

class StructureBasedModel(Model):

    def __init__(self, topology, bead_repr=None):
        Model.__init__(self, topology, bead_repr=bead_repr)
        self.Hamiltonian = ptl.StructureBasedHamiltonian()
         
    def set_reference(self, traj):
        self.ref_traj_aa = traj[0]
        self.ref_traj = self.mapping.map_traj(traj[0])

    def add_sbm_potentials(self):
        self.Hamiltonian.add_sbm_potentials(self) 
