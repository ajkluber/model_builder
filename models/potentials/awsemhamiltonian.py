import numpy as np

import mdtraj as md

import util

import awsem

from hamiltonian import Hamiltonian

class StructureBasedHamiltonian(Hamiltonian):

    def __init__(self):
        Hamiltonian.__init__(self)
        self.gamma_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
                               'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                               'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
                               'SER', 'THR', 'TRP', 'TYR', 'VAL']

    def _source_parameters(self, param_path, fix_backbone_path):
        self._param_path = param_path 
        self._source_gammas()
        self._source_fix_backbone_coeff(fix_backbone_coeff_path)

    def _source_fix_backbone_coeff(self, fix_backbone_coeff_path):
        pass

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

    def calculate_local_density(traj):
        #TODO
        # all contact pairs
        contact_pairs = []
        for i in range(traj.n_residues):
            if traj.top.residue(i).name == "GLY":
                idx1 = traj.top.select("resid {} and name CA".format(i))[0]
            else:
                idx1 = traj.top.select("resid {} and name CB".format(i))[0]
            for j in range(traj.n_residues):
                if i != j:
                    if traj.top.residue(j).name == "GLY":
                        idx2 = traj.top.select("resid {} and name CA".format(j))[0]
                    else:
                        idx2 = traj.top.select("resid {} and name CB".format(j))[0]
                    contact_pairs.append([idx1, idx2])
        contact_pairs = np.array(contact_pairs)

        contact_dists = md.compute_distances(traj, contact_pairs)[0]
        contacts = awsem.theta_I(contact_dists)
        contact_matrix = contacts.reshape((traj.n_residues, traj.n_residues - 1))
        local_density = np.sum(contact_matrix, axis=1)
        return local_density

    def calculate_energy(self):
        pass

