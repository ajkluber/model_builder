import numpy as np

from model_builder.models.structure import contacts as cts


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

class Potentials(object):
    """Mode Hamiltonian"""
    pair_potentials = {2:LJ1210Potential}

    def __init__(self):
        self._bonds = []
        self._angles = []
        self._dihedrals = []
        self._pairs = []

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
            print "warning"

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

############################################################################
# Pair potentials
############################################################################

class PairPotential(object):
     
    def __init__(self, atmi, atmj):
        self.atmi = atmi
        self.atmj = atmj

    def describe(self):
        """interaction description"""
        return "%s:%12s%12s" % (self.prefix_label, self.atmi, self.atmj)

    def __hash__(self):
        hash_value = hash(self.prefix_label)
        hash_value ^= hash(self.atmi.index)
        hash_value ^= hash(self.atmj.index)
        return hash_value

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

class LJPotential(PairPotential):
    
    def __init__(self, atmi, atmj, eps, r0):
        PairPotential.__init__(self, atmi, atmj)
        self.eps = eps
        self.r0 = r0

    def __hash__(self):
        hash_value = PairPotential.__hash__(self)
        hash_value ^= hash(self.eps)
        hash_value ^= hash(self.r0)
        return hash_value

class LJ12Potential(LJPotential):

    def __init__(self, atmi, atmj, eps, r0):
        LJPotential.__init__(self, atmi, atmj, eps, r0) 
        self.prefix_label = "LJ12"

    def V(self, r):
        return self.eps*self.dVdeps(r)

    def dVdeps(self, r):
        x = self.r0/r
        return x**12

    def dVdr(self, r):
        return self.eps*self.d2Vdrdeps(r)

    def d2Vdrdeps(self, r):
        x = self.r0/r
        return (-60./self.r0)*((x**13) - (x**11))

class LJ1210Potential(LJPotential):
    
    def __init__(self, atmi, atmj, r0, eps):
        LJPotential.__init__(self, atmi, atmj, r0, eps)
        self.prefix_label = "LJ1210"

    def V(self, r):
        return self.eps*self.dVdeps(r)

    def dVdeps(self, r):
        x = self.r0/r
        return 5.*(x**12) - 6.*(x**10)

    def dVdr(self, r):
        return self.eps*self.d2Vdrdeps(r)

    def d2Vdrdeps(self, r):
        x = self.r0/r
        return (-60./self.r0)*((x**13) - (x**11))

class TanhRepPotential(PairPotential):
    
    def __init__(self, atmi, atmj, eps, r0, width):
        PairPotential.__init__(self, atmi, atmj)
        self.prefix_label = "TANHREP"
        self.eps = eps
        self.r0 = r0
        self.width = width

    def V(self,r):
        return self.eps*dVdeps(r) 

    def dVdeps(self, r):
        alpha = 1./self.width
        r0prime = self.r0 + self.width
        return 0.5*(np.tanh(-alpha*(r - r0prime)) + 1.)

    def dVdr(self, r):
        return self.eps*d2Vdrdeps(r) 

    def d2Vdrdeps(self, r):
        alpha = 1./self.width      
        r0prime = self.r0 + self.width
        return -0.5*alpha*(1. - (np.tanh(-alpha*(r - r0prime)))**2)

    def __hash__(self):
        hash_value = PairPotential.__hash__(self)
        hash_value ^= hash(self.eps)
        hash_value ^= hash(self.r0)
        hash_value ^= hash(self.width)
        return hash_value

if __name__ == "__main__":
    pass
