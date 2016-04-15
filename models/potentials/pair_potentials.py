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
    
    def _list_hash(self):
        listhash = [hash(self.prefix_label)]
        listhash.append(hash(self.atmi))
        listhash.append(hash(self.atmj))
        
        return listhash
        
    def __eq__(self, other):
        test = True
        for i,j in zip(self._list_hash(), other._list_hash()):
            test = test and (i==j)
            
        return test

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
        return (-12./self.r0)*(x**13)

class LJ126Potential(LJPotential):

    def __init__(self, atmi, atmj, eps, r0):
        LJPotential.__init__(self, atmi, atmj, eps, r0) 
        self.prefix_label = "LJ126"

    def V(self, r):
        return self.eps*self.dVdeps(r)

    def dVdeps(self, r):
        x = self.r0/r
        return 4*((x**12) - (x**6))

    def dVdr(self, r):
        return self.eps*self.d2Vdrdeps(r)

    def d2Vdrdeps(self, r):
        x = self.r0/r
        return (-24/self.r0)*(2.*(x**13) - (x**7))

class LJ1210Potential(LJPotential):
    
    def __init__(self, atmi, atmj, eps, r0):
        LJPotential.__init__(self, atmi, atmj, eps, r0)
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

class LJ1210RepPotential(LJPotential):
    
    def __init__(self, atmi, atmj, eps, r0):
        LJPotential.__init__(self, atmi, atmj, eps, r0)
        self.prefix_label = "LJ1210REP"

    def V(self, r):
        return self.eps*self.dVdeps(r)

    def dVdeps(self, r):
        x = self.r0/r
        V = np.zeros(x.shape,float)
        V[x > 1] = 5.*(x[x > 1]**12) - 6.*(x[x > 1]**10) + 2.
        V[x <= 1] = -5.*(x[x <= 1]**12) - 6.*(x[x <= 1]**10)
        return V

    def dVdr(self, r):
        return self.eps*self.d2Vdrdeps(r)

    def d2Vdrdeps(self, r):
        x = self.r0/r
        V = np.zeros(x.shape,float)
        V[x > 1] = (-60./r0)*(x[x > 1]**13 - x[x > 1]**11)
        V[x <= 1] = (60./r0)*(x[x <= 1]**13 - x[x <= 1]**11)
        return V

class TanhRepPotential(PairPotential):
    
    def __init__(self, atmi, atmj, eps, r0, width):
        PairPotential.__init__(self, atmi, atmj)
        self.prefix_label = "TANHREP"
        self.eps = eps
        self.r0 = r0
        self.width = width

    def V(self,r):
        return self.eps*self.dVdeps(r) 

    def dVdeps(self, r):
        alpha = 1./self.width
        r0prime = self.r0 + self.width
        return 0.5*(np.tanh(-alpha*(r - r0prime)) + 1.)

    def dVdr(self, r):
        return self.eps*self.d2Vdrdeps(r) 

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

class LJ12TanhRepPotential(PairPotential):
    
    def __init__(self, atmi, atmj, eps, rNC, r0, width):
        PairPotential.__init__(self, atmi, atmj)
        self.prefix_label = "LJ12TANHREP"
        self.eps = eps
        self.rNC = rNC
        self.r0 = r0
        self.width = width

    def V(self,r):
        return self.eps*self.dVdeps(r) + (self.rNC/r)**12

    def dVdeps(self, r):
        alpha = 1./self.width
        r0prime = self.r0 + self.width
        return 0.5*(np.tanh(-alpha*(r - r0prime)) + 1.)

    def dVdr(self, r):
        return self.eps*self.d2Vdrdeps(r) - (12./self.rNC)*((self.rNC/r)**13)

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

    def dVdeps(self, r):
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

class LJ12GaussianPotential(PairPotential):

    def __init__(self, atmi, atmj, eps, rNC, r0, width):
        PairPotential.__init__(self, atmi, atmj)
        self.prefix_label = "LJ12GAUSSIAN"
        self.eps = eps
        self.rNC = rNC
        self.r0 = r0
        self.width = width
        self.gaussian = GaussianPotential(atmi, atmj, eps, r0, width)
        self.lj12 = LJ12Potential(atmi, atmj, 1.0, rNC)

    def V(self, r):
        return (1. + self.lj12.V(r))*(1. + self.gaussian.V(r)) - 1.

    def dVdr(self, r):
        first = (1. + self.lj12.V(r))*self.gaussian.dVdr(r)
        second = (1. + self.gaussian.V(r))*self.lj12.dVdr(r)
        return first + second

    def dVdeps(self, r):
        return self.gaussian.dVdeps(r)
        
    def d2Vdrdeps(self, r):
        return self.gaussian.d2Vdrdeps(self, r)

    def __hash__(self):
        hash_value = PairPotential.__hash__(self)
        hash_value ^= hash(self.eps)
        hash_value ^= hash(self.rNC)
        hash_value ^= hash(self.r0)
        hash_value ^= hash(self.width)
        return hash_value

class FlatBottomWell(PairPotential):

    def __init__(self, atmi, atmj, kb, rNC, r0):
        PairPotential.__init__(self, atmi, atmj)
        self.prefix_label = "FLATWELL"
        self.kb = kb
        self.rNC = rNC
        self.r0 = r0

    def V(self, r):
        V = np.empty(r.shape[0])
        V[r < self.r0] = (rNC/r[r < self.r0])**12
        V[r >= self.r0] = 0.5*(r[r >= self.r0] - self.r0)**2
        return V

    def V(r, kb, rNC, r0):
        V = np.empty(r.shape[0])
        V[r < r0] = (rNC/r[r < r0])**12
        V[r >= r0] = 0.5*kb*((r[r >= r0] - r0)**2)
        return V

PAIR_POTENTIALS = {"LJ1210":LJ1210Potential,
                "GAUSSIAN":GaussianPotential,
                "LJ12GAUSSIAN":LJ12GaussianPotential,
                "TANHREP":TanhRepPotential,
                "LJ12TANHREP":LJ12TanhRepPotential,
                "FLATWELL":FlatBottomWell}
