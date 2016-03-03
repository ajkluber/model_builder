import numpy as np

# TODO: Add derivatives of potential w.r.t. all parameters

############################################################################
# Pair potentials
############################################################################

class PairPotential(object):
     
    def __init__(self, atmi, atmj):
        self.atmi = atmi
        self.atmj = atmj

    def describe(self):
        """interaction description"""
        return "{}:{:>12}{:>12}".format(self.prefix_label, self.atmi, self.atmj)

    def __hash__(self):
        hash_value = hash(self.prefix_label)
        hash_value ^= hash(self.atmi)
        hash_value ^= hash(self.atmj)
        return hash_value

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __repr__(self):
        return "<PairPotential at 0x{}x>".format(id(self))

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

class GaussianPotential(PairPotential):

    def __init__(self, atmi, atmj, eps, r0, width):
        PairPotential.__init__(self, atmi, atmj)
        self.prefix_label = "GAUSSIAN"
        self.eps = eps
        self.r0 = r0
        self.width = width

    def V(self, r):
        return self.eps*self.dVdeps(r)

    def dVdeps(self, r)
        return -np.exp(-((r - self.r0)**2)/(2.*(self.width**2)))

    def dVdr(self, r):
        return self.eps*self.d2Vdrdeps(r)

    def d2Vdrdeps(self, r):
        return ((r - self.r0)/(self.width**2))*np.exp(-((r - self.r0)**2)/(2.*(self.width**2)))

    def __hash__(self):
        hash_value = PairPotential.__hash__(self)
        hash_value ^= hash(self.eps)
        hash_value ^= hash(self.r0)
        hash_value ^= hash(self.width)
        return hash_value

class LJ12GaussianPotential(PairPotential,LJ12Potential,GaussianPotential):

    def __init__(self, atmi, atmj, eps, rNC, r0, width):
        PairPotential.__init__(self, atmi, atmj)
        self.prefix_label = "LJ12GAUSSIAN"
        self.eps = eps
        self.rNC = rNC
        self.r0 = r0
        self.width = width

    def V(self, r):
        return LJ12Potential.V(self, r)*(1. + GaussianPotential.dVdeps(self, r))

    def dVdr(self, r):
        return LJ12Potential.dVdr(self, r)*(1. + GaussianPotential.dVdeps(self, r)) +\
               LJ12Potential.dVdr(self, r)*GaussianPotential.d2Vdrdeps(self, r)

    def dVdeps(self, r):
        return LJ12Potential.dVdeps(self, r)*(1. + GaussianPotential.dVdeps(self, r))

    def d2Vdrdeps(self, r):
        return LJ12Potential.dVdr(self, r)*(1. + GaussianPotential.dVdeps(self, r)) +\
               LJ12Potential.dVdr(self, r)*GaussianPotential.d2Vdrdeps(self, r)

    def __hash__(self):
        hash_value = PairPotential.__hash__(self)
        hash_value ^= hash(self.eps)
        hash_value ^= hash(self.rNC)
        hash_value ^= hash(self.r0)
        hash_value ^= hash(self.width)
        return hash_value

PAIR_POTENTIALS = {"LJ1210":LJ1210Potential,
                "GAUSSIAN":GaussianPotential,
                "LJ12GAUSSIAN":LJ12GaussianPotential,
                "TANHREP":TanhRepPotential}
