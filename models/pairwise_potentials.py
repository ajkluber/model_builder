""" Utilities for nonbonded/pairwise potential terms.


"""

import numpy as np

## To Do:
## - Dictionary of pairwise potential types

def wrap_pairwise(old_pairwise,*args):
    """ Wraps pairwise function so only r needs to be passed

    Positional arguments need to be passed correctly.
    """
    def new_pairwise(r):
        return old_pairwise(r,*args)
    return new_pairwise

def LJ12(r,sigma):
    x = r/sigma
    V = x**12
    return V

def LJ12_deriv(r,sigma):
    x = r/sigma
    V = (-12./sigma)*(x**13)
    return V

def LJ1210(r,sigma):
    x = r/sigma
    V = 5*(x**12) - 6*(x**10)
    return V

def LJ1210_deriv(r,sigma):
    x = r/sigma
    V = (-60./sigma)*((x**13) - (x**11))
    return V

def LJ1210rep(r,sigma):
    x = r/sigma
    V = np.zeros(x.shape,float)
    V[x > 1] = (5*(x[x > 1]**12) - 6*(x[x > 1]**10) + 2)
    V[x <= 1] = -(5*(x[x <= 1]**12) - 6*(x[x <= 1]**10))
    return V

def LJ1210rep_deriv(r,sigma):
    x = r/sigma
    V = np.zeros(x.shape,float)
    V[x > 1] = -60*(1./sigma)*(x[x > 1]**13 - x[x > 1]**11)
    V[x <= 1] = 60*(1./sigma)*(x[x <= 1]**13 - x[x <= 1]**11)
    return V

def Gaussian(r,r0,width):
    V = -np.exp(-((r - r0)**2)/(2.*(width**2)))
    return V

def Gaussian_deriv(r,r0,width):
    V = ((r - r0)/(width**2))*np.exp(-((r - r0)**2)/(2.*(width**2)))
    return V
