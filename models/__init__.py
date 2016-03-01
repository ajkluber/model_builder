""" Submodule with classes of coarse-grain models.

Description:

    A module that prepares the Model object which contains the functions and
parameters needed to prepare input files for coarse-grain protein simulations
in Gromacs. Each model is a seperate class. 


Classes:

References:

"""
import numpy as np

import CoarseGrainedModel
import bonded_potentials
import pairwise_potentials
import pdb_parser
import prep_Gaussian_pairwise


import structure
import potentials
import output

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

    def map_traj(self, traj):
        self.structure_mapping.map_traj(traj)

