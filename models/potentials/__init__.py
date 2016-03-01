import numpy as np

from model_builder.models.structure import contacts as cts
from model_builder.models.potentials.pair_potentials import *


#POTENTIALS = {1:LJ12Potential,
#                 2:LJ1210LJ12Potential,
#                 3:LJ1210repLJ12Potential,
#                 4:GaussianLJ12Potential,
#                 5:Cheng_repLJ12Potential,
#                 6:LJ126LJ12Potential,
#                 7:LJ126repLJ12Potential,
#                 8:compound_LJ12_GaussianLJ12Potential,
#                 9:compound_LJ12_double_GaussianLJ12Potential,
#                 10:compound_double_GaussianLJ12Potential,
#                 11:FRET_EfficiencyLJ12Potential}[code]


def _hash_numpy_array(x):
    hash_value = hash(x.shape)
    hash_value ^= hash(x.strides)
    hash_value ^= hash(x.data.tobytes())
    return hash_value

class Hamiltonian(object):
    """Mode Hamiltonian"""
    pair_potentials = \
        {2:LJ1210Potential
        5:TanhRepPotential
        }

    def __init__(self):
        self._bonds = []
        self._angles = []
        self._dihedrals = []
        self._pairs = []
        #self._special = []  # Add?

    @property
    def n_bonds(self):
        return len(self._bonds)

    @property
    def n_angles(self):
        return len(self._angles)

    @property
    def n_dihedrals(self):
        return len(self._dihedrals)

    @property
    def n_pairs(self):
        return len(self._pairs)

    @property
    def bonds(self):
        for bond in self._bonds:
            yield bond

    @property
    def angles(self):
        for angle in self._angles:
            yield angle

    @property
    def dihedrals(self):
        for dihedral in self._dihedrals:
            yield dihedral

    @property
    def pairs(self):
        for pair in self._pairs:
            yield pair

    @property
    def potentials(self):
        pots = self._bonds + self._angles + self._dihedrals + self._pairs
        for pot in pots:
            yield pot 

    def describe(self):
        labels = []
        for pot in self.pairs:
            labels.append(pot.describe())
        return labels

    def add_pair(self, code, atm1, atm2, *args):
        p = pair_potentials[code](atm1, atm2, *args)
        if p not in self._pairs:
            self._pairs.append(p)
        else:
            print "warning: pair already has this interaction {}. skipping.".format(p.describe())

    def add_sbm_contacts(self, Model):
        residue_contacts = cts.residue_contacts(Model.ref_traj)
        atm_pairs = Model.structure_mapping.residue_to_atom_contacts(residue_contacts)

        code = 2
        eps = 1.
        xyz = Model.ref_traj.xyz[0]
        self.pairV = []
        for atm1, atm2 in atm_pairs:
            r0 = np.linalg.norm(xyz[atm1.index,:] - xyz[atm2.index,:])
            self.add_pair(code, atm1, atm2, eps, r0)

    def define_contact_group(self, label, pairs):
        # Use this to define a group of contacts by a label.
        # The label can later be used to get their energy. 
        # e.g. 'native' [[1,10],[2,10]]
        pass 

    def calc_contact_group_energy(self, label, traj):
        # Calculate the energy of a group of contacts
        pass

    def select_parameters(self):
        # Identify parameters by:
        #   - parameter type: eps, r0, 
        #   - interaction type: bonds, angles, etc.
        pass


if __name__ == "__main__":
    pass
