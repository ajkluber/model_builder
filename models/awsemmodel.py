import numpy as np

import potentials 

from model import Model            


class AwsemModel(Model):

    def __init__(self, topology):
        Model.__init__(self, topology, bead_repr="AWSEM")
        self.Hamiltonian = potentials.AwsemHamiltonian(self.mapping.top)
        self.Hamiltonian._set_charged_residues(self.mapping._charged_residues)

    def source_parameters(self, param_path, parameterize=True): 
        # read in parameters from source directory
        self.Hamiltonian._source_parameters(param_path)

