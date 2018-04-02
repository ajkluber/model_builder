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

    def set_epsilon(self, value):
        self.eps = value

    def get_V_epsilons(self, r):
        """ Returns function V(epsilons)

        Default function for returning the Potential Energy as a
        function of epsilons. Since the majority of pairwise functions
        are scaled linearly with epsilon, this is a good default choice. This method can and should be overrided when necessary. See

        Parameters
        ----------
        r : array(float)
            Distance for evaluating each pairwise potential function.

        Returns
        -------
        func : method
            Function that computse the potential energy as a function
            of epsilon
        """

        constants_list = self.dVdeps(r)
        def func(epsilon):
            return constants_list * epsilon

        return func

    def get_dV_depsilons(self, r):
        """ Returns function dV(epsilons)/depsilons

        Parameters
        ----------
        r : array(float)
            Distance for evaluating each pairwise potential function.

        Returns
        -------
        func : method
            Function that computse the derivative of the potential
            energy with respect to epsilon, as a function of epsilon.


        """

        constants_list = self.dVdeps(r)
        def func(epsilon):
            return constants_list

        return func

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
        self.other_params = [r0]

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
        self.other_params = [r0, width]

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

class LJ12TanhRepPotential(PairPotential):

    def __init__(self, atmi, atmj, eps, rNC, r0, width):
        PairPotential.__init__(self, atmi, atmj)
        self.prefix_label = "LJ12TANHREP"
        self.eps = eps
        self.rNC = rNC
        self.r0 = r0
        self.width = width
        self.other_params = [rNC, r0, width]

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

class GaussianPotential(PairPotential):

    def __init__(self, atmi, atmj, eps, r0, width):
        PairPotential.__init__(self, atmi, atmj)
        self.prefix_label = "GAUSSIAN"
        self.eps = eps
        self.r0 = r0
        self.width = width
        self.other_params = [r0, width]

    def V(self, r):
        return self.eps*self.dVdeps(r)

    def dVdeps(self, r):
        return -np.exp(-((r - self.r0)**2)/(2.*(self.width**2)))

    def dVdr(self, r):
        return self.eps*self.d2Vdrdeps(r)

    def d2Vdrdeps(self, r):
        return ((r - self.r0)/(self.width**2))*np.exp(-((r - self.r0)**2)/(2.*(self.width**2)))

class LJ12GaussianPotential(PairPotential):

    def __init__(self, atmi, atmj, eps, rNC, r0, width):
        PairPotential.__init__(self, atmi, atmj)
        self.prefix_label = "LJ12GAUSSIAN"
        self.eps = eps
        self.rNC = rNC
        self.r0 = r0
        self.width = width
        self.gaussian = GaussianPotential(atmi, atmj, 1.0, r0, width)
        self.lj12 = LJ12Potential(atmi, atmj, 1.0, rNC)
        self.other_params = [rNC, r0, width]

    def V(self, r):
        return (self.eps * self.gaussian.V(r)) + self.lj12.V(r) + (self.lj12.V(r) * self.gaussian.V(r))

    def dVdr(self, r):
        return (self.eps * self.gaussian.dVdr(r)) + self.lj12.dVdr(r) + (self.lj12.dVdr(r) * self.gaussian.V(r))+ (self.lj12.V(r) * self.gaussian.dVdr(r))

    def dVdeps(self, r):
        return self.gaussian.dVdeps(r)

    def d2Vdrdeps(self, r):
        return self.gaussian.d2Vdrdeps(r)

    def set_epsilon(self, value):
        self.eps = value

class LJ12GaussTanhSwitching(PairPotential):
    """ LJ12 Potential with Gaussian attractive and tanh repulsive"""
    def __init__(self, atmi, atmj, eps, rNC, r0, width):
        PairPotential.__init__(self, atmi, atmj)
        self.prefix_label = "LJ12GAUSSIANTANH"
        self.eps = eps
        self.rNC = rNC
        self.r0 = r0
        self.width = width
        self.attractive = LJ12GaussianPotential(atmi, atmj, np.abs(eps), rNC, r0, width)
        self.repulsive = LJ12TanhRepPotential(atmi, atmj, np.abs(eps), rNC, r0, width)
        self.determine_current()
        self.other_params = [rNC, r0, width]

    def V(self, r):
        return self.current.V(r)

    def dVdr(self, r):
        return self.current.dVdr(r)

    def dVdeps(self, r):
        return self.current.dVdeps(r)

    def d2Vdrdeps(self, r):
        return self.current.d2Vdrdeps(self, r)

    def determine_current(self):
        """ If eps > 0, return attractive, otherwise return repulsive"""
        if self.eps < 0:
            self.current = self.repulsive
        else:
            self.current = self.attractive

    def set_epsilon(self, value):
        self.eps = value
        self.attractive.set_epsilon(np.abs(value))
        self.repulsive.eps = np.abs(value)
        self.determine_current()

    def get_V_epsilons(self, r):
        constants_list_att = self.attractive.dVdeps(r)
        constants_list_rep = self.repulsive.dVdeps(r)
        def func(epsilon):
            if epsilon < 0:
                return constants_list_rep * epsilon
            else:
                return constants_list_att * epsilon

        return func

    def get_dV_depsilons(self, r):
        constants_list_att = self.attractive.dVdeps(r)
        constants_list_rep = self.repulsive.dVdeps(r)
        constants_list_average = (constants_list_att + constants_list_rep) / 2.
        def func(epsilon):
            if epsilon < 0:
                return constants_list_rep
            elif epsilon == 0:
                return constants_list_average
            else:
                return constants_list_att

        return func



class FlatBottomWell(PairPotential):

    def __init__(self, atmi, atmj, kb, rNC, r0):
        PairPotential.__init__(self, atmi, atmj)
        self.prefix_label = "FLATWELL"
        self.kb = kb
        self.rNC = rNC
        self.r0 = r0
        self.other_params = [kb, rNC, r0]

    def V(self, r):
        V = np.zeros(r.shape[0])
        V[r < self.r0] = (self.rNC/r[r < self.r0])**12
        V[r >= self.r0] = 0.5*self.kb*((r[r >= self.r0] - self.r0)**2)
        return V

    def dVdr(self, r):
        dVdr = np.zeros(r.shape[0])
        dVdr[r < self.r0] = -(12./self.rNC)*((self.rNC/r[r < self.r0])**13)
        dVdr[r >= self.r0] = self.kb*(r[r >= self.r0] - self.r0)
        return dVdr

class CustomPairPotential(PairPotential):

    def __init__(self, atmi, atmj, func, *args):
        PairPotential.__init__(self, atmi, atmj)
        self.prefix_label = "CUSTOM"
        self.args = args
        self.func = func

    def V(self, r):
        return self.func(r, *self.args)

    def dVdr(self, r):
        return np.gradient(self.func(r, *self.args), r[1] - r[0])

PAIR_POTENTIALS = {"LJ1210":LJ1210Potential,
                "GAUSSIAN":GaussianPotential,
                "LJ12GAUSSIAN":LJ12GaussianPotential,
                "TANHREP":TanhRepPotential,
                "LJ12TANHREP":LJ12TanhRepPotential,
                "CUSTOM":CustomPairPotential,
                "FLATWELL":FlatBottomWell,
                "LJ12GAUSSIANTANH":LJ12GaussTanhSwitching}
