import warnings
import numpy as np

import mdtraj as md

from model_builder.models.structure import contacts as cts
#from model_builder.models.potentials.pair_potentials import *
#from model_builder.models.potentials.bonded_potentials import *
import pair_potentials 
import bonded_potentials


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

def interaction_exists_warning(pot):
    warnings.warn("Interaction already exists! skipping: {}".format(pot.describe()))

def default_sbm_parameters_warning():
    warning.warn("Using default SBM parameters")

def default_sbm_potentials_warning():
    warning.warn("Using default SBM parameters")

def missing_reference_warning():
    warning.warn("Need to set reference structure model.set_reference()")

class Hamiltonian(object):
    """Mode Hamiltonian"""

    def __init__(self):
        self._bonds = []
        self._angles = []
        self._dihedrals = []
        self._pairs = []
        self._default_parameters = {}
        self._default_potentials = {}

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
    
    def _add_bond(self, code, atm1, atm2, *args):
        b = bonded_potentials.BOND_POTENTIALS[code](atm1, atm2, *args)
        if b not in self._bonds:
            self._bonds.append(b)
        else:
            interaction_exists_warning(b)

    def _add_angle(self, code, atm1, atm2, atm3, *args):
        ang = bonded_potentials.ANGLE_POTENTIALS[code](atm1, atm2, atm3, *args)
        if ang not in self._angles:
            self._angles.append(ang)
        else:
            interaction_exists_warning(ang)

    def _add_dihedral(self, code, atm1, atm2, atm3, atm4, *args):
        dih = bonded_potentials.DIHEDRAL_POTENTIALS[code](atm1, atm2, atm3, atm4, *args)
        if dih not in self._dihedrals:
            self._dihedrals.append(dih)
        else:
            interaction_exists_warning(dih)

    def _add_pair(self, code, atm1, atm2, *args):
        p = pair_potentials.PAIR_POTENTIALS[code](atm1, atm2, *args)
        if p not in self._pairs:
            self._pairs.append(p)
        else:
            interaction_exists_warning(p)

    def _add_bonds(self, bond_params):
        for b in bond_params:
            self._add_bond(b[0], b[1], b[2], b[3:])

    def _add_angles(self, angle_params):
        for a in angle_params:
            self._add_angle(a[0], a[1], a[2], a[3], a[4:])

    def _add_dihedrals(self, dihedral_params):
        for d in dihedral_params:
            self._add_dihedral(d[0], d[1], d[2], d[3], d[4], d[5:])

    def _add_pairs(self, pair_params):
        for p in pair_params:
            self._add_pair(p[0], p[1], p[2], p[3:])

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

    def set_parameters(self):
        # Some way to set parameters?
        pass

class StructureBasedHamiltonian(Hamiltonian):

    def __init__(self):
        Hamiltonian.__init__(self)

    def add_sbm_potentials(self, Model):

        if not hasattr(Model,"ref_traj"):
            raise AttributeError("Need to set reference structure model.set_reference()")

        # How should user defined parameters be set? 
        if self._default_parameters == []:
            default_sbm_parameters_warning()
            self._use_sbm_default_parameters() 

        if self._default_potentials == []:
            default_sbm_potentials_warning()
            self._use_sbm_default_potentials() 

        Model.mapping._assign_sbm_angles()
        Model.mapping._assign_sbm_dihedrals()
        Model.mapping._assign_sbm_contacts()

        self._add_sbm_bonds(Model)
        self._add_sbm_angles(Model)
        self._add_sbm_dihedrals(Model)
        self._add_sbm_contacts(Model)

    def _add_sbm_bonds(self, Model):
        top = Model.mapping.top

        if hasattr(Model,"ref_traj"):
            kb = self._default_parameters["kb"]
            code = self._default_potentials["bond"]
            xyz = Model.ref_traj.xyz[0]
            for atm1, atm2 in top.bonds:
                r0 = np.linalg.norm(xyz[atm1.index,:] - xyz[atm2.index,:])
                self._add_bond(code, atm1, atm2, kb, r0)
        else:
            missing_reference_warning()

    def _add_sbm_angles(self, Model):
        structure = Model.mapping

        if hasattr(Model,"ref_traj"):
            ka = self._default_parameters["ka"]
            code = self._default_potentials["angles"]
            for atm1, atm2, atm3 in structure._angles:
                theta0 = md.compute_anlges(Model.ref_traj, \
                        atm1.index, atm2.index, atm3.index)[0][0]
                self._add_angle(code, atm1, atm2, atm3, ka, theta0)
        else:
            missing_reference_warning()

    def _add_sbm_dihedrals(self, Model):

        structure = Model.mapping

        if hasattr(Model,"ref_traj"):
            kd = self._default_parameters["kd"]
            code = self._default_potentials["dihedrals"]
            for atm1, atm2, atm3 in structure._dihedrals:
                phi0 = 180. + (180./np.pi)*md.compute_dihedral(Model.ref_traj, \
                            atm1.index, atm2.index, atm3.index, atm4.index)[0][0]
                self._add_angle(code, atm1, atm2, atm3, atm4, kd, phi0)
        else:
            missing_reference_warning()

    def _add_sbm_contacts(self, Model):
        """Add structure-based model contacts"""
        # TODO: Allow for different contact code's (e.g. LJ1210, Gaussian, etc.)
    
        if self._default_parameters == []:
            default_sbm_parameters_warning()
            self._use_sbm_default_parameters() 

        residue_contacts = cts.residue_contacts(Model.ref_traj)
        atm_pairs = Model.mapping._residue_to_atom_contacts(residue_contacts)

        code = 2    # LJ1210 for now
        if hasattr(Model,"ref_traj"):
            eps = self._default_parameters["eps"]
            xyz = Model.ref_traj.xyz[0]
            for atm1, atm2 in atm_pairs:
                r0 = np.linalg.norm(xyz[atm1.index,:] - xyz[atm2.index,:])
                self._add_pair(code, atm1, atm2, eps, r0)
        else:
            missing_reference_warning()

    def _use_sbm_default_parameters(self):
        self._default_parameters = {"kb":20000., "ka":40.,
                                    "kd":1., "eps":1}

    def _use_sbm_default_potentials(self):
        self._default_potentials = {"bond":"HARMONIC_BOND",
                                    "angle":"HARMONIC_ANGLE",
                                    "dihedral":"COSINE_DIHEDRAL"}
        

if __name__ == "__main__":
    pass
