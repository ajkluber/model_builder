
#import model_builder.models.interactions.pair_potentials as pair_potentials



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

#__all__ = ["Potentials", "PairPotential"]

class Potentials(object):

    def __init__(self):
        self._pairs = []
        self._bonds = []
        self._angles = []
        self._dihedrals = []

    def add_pair(self, code, atm1, atm2, *args):
        self._pairs.append(POTENTIALS[code](atm1, atm2, *args))

    @property
    def n_pairs(self):
        return len(self._pairs)

    @property
    def n_bonds(self):
        return len(self._bonds)

    @property
    def n_angles(self):
        return len(self._angles)

    @property
    def n_dihedrals(self):
        return len(self._dihedrals)

    def describe(self):
        pass

class PairPotential(object):
    
    def __init__(self, atmi, atmj):
        self.atmi = atmi
        self.atmj = atmj

    def __repr__(self):
        """interaction description"""
        return "%s:%12s%12s" % (self.prefix_label, self.atmi, self.atmj)

class LJ1210Potential(PairPotential):
    
    def __init__(self, atmi, atmj, r0, eps):
        PairPotential.__init__(self, atmi, atmj)
        #self.code = 2
        self.prefix_label = "LJ1210"
        self.r0 = r0
        self.eps = eps

    def V(self, r):
        return self.eps*self.dim_V(r)

    def dim_V(self, r):
        x = self.r0/r
        V = 5.*(x**12) - 6.*(x**10)
        return V

    def dVdr(self, r):
        return self.eps*self.dim_dVdr(r)

    def dim_dVdr(self, r):
        x = self.r0/r
        V = (-60./self.r0)*((x**13) - (x**11))
        return V
    
POTENTIALS = {2:LJ1210Potential}

if __name__ == "__main__":
    pass
