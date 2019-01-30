from __future__ import absolute_import
import numpy as np

import model_builder.models.mappings as mappings
import model_builder.models.potentials as potentials

class Model(object):
    """Model class """
    def __init__(self, topology, bead_repr="CA"):
        self.mapping = mappings.assign_mapping(bead_repr, topology)
        self.Hamiltonian = potentials.Hamiltonian()

    def describe(self):
        #TODO: What is the best description?
        pass

    def map_traj(self, traj):
        return self.mapping.map_traj(traj)

    def save_starting_conf(self, saveas="conf.gro"):
        self.starting_traj.save(saveas)

    def set_starting_conf(self, traj):
        self.starting_traj = traj

    def add_pairs(self, pairs):
        self.mapping._add_pairs(pairs)

    def assign_fitted_epsilons(self, params_to_fit_indices):
        self.params_to_fit_indices = params_to_fit_indices
        self.fitted_epsilons = []
        self.fitted_function_types = []
        for i in params_to_fit_indices:
            self.fitted_epsilons.append(self.Hamiltonian._epsilons[i])
            self.fitted_function_types.append(self.Hamiltonian._pair_function_type_labels[i])

    def output_epsilons(self):
        params = self.Hamiltonian._epsilons
        for idx, i in enumerate(self.params_to_fit_indices):
            params[i] = self.fitted_epsilons[idx]

        f = open("params", "w")
        f.write("# Fitted Parameters\n")
        for param in params:
            f.write("%f\n"%param)
        f.close()
