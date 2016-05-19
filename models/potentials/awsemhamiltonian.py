import numpy as np

import mdtraj as md

import util

from hamiltonian import Hamiltonian

class StructureBasedHamiltonian(Hamiltonian):

    def __init__(self):
        Hamiltonian.__init__(self)
        self.gamma_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
                               'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                               'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
                               'SER', 'THR', 'TRP', 'TYR', 'VAL']

    def _source_parameters(self, param_path):

        self._param_path = param_path 
        self._source_gammas()


    def _source_gammas(self):
        with open("{}/gamma.dat".format(self.param_path), "r") as fin:
            gamma_direct = np.zeros((21,20))
            for i in range(21):
                for j in range(i, 20):
                    line = fin.readline()
                    gamma_direct[i, j] = float(line.split()[0])

            line = fin.readline()
            gamma_water = np.zeros((21,20))
            gamma_protein = np.zeros((21,20))
            for i in range(21):
                for j in range(20):
                    line = fin.readline()
                    gamma_water[i, j] = float(line.split()[0])
                    gamma_protein[i, j] = float(line.split()[1])

        self.gamma_direct = gamma_direct 
        self.gamma_water = gamma_water 
        self.gamma_protein = gamma_protein 

    def calculate_energy(self):
        pass
