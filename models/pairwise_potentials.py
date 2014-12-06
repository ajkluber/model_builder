""" Utilities for nonbonded/pairwise potential terms.


"""

import numpy as np

## To Do:
## - Dictionary of pairwise potential types. extensible

def get_pair_potential(code):
    """ Returns pairwise potential function"""
    potential = {1:LJ12,2:LJ1210,3:LJ1210rep,4:Gaussian}[code]
    return potential

def get_pair_potential_deriv(code):
    """ Returns derivative of pairwise potential function"""
    potential = {1:LJ12_deriv,2:LJ1210_deriv,3:LJ1210rep_deriv,4:Gaussian_deriv}[code]
    return potential

def wrap_pairwise(old_pairwise,*args):
    """ Wraps pairwise function so only r needs to be passed 

    Inputs:
        A function that takes one distance (r) and a variable number of
    positional arguments depends on (such as equilibrium distance or 
    width of well).

    Returns:
        A function that only depends on the distance, customized to fixed values
    of the remaining positional parameters.
    """
    def new_pairwise(r):
        return old_pairwise(r,*args)
    return new_pairwise

def wrap_sum_of_pairwise(*args):
    """ Returns function that is sum of args functions. 

    Inputs:
        A list of functions that each depend on a distance (r).
    Returns:
        A function that evaluates the sum of the inputted functions.
    The returned function is also depends on the distance (r).
    NOT USED
    """
    def sum_pairwise(r):
        return sum(np.array([ args[i](r) for i in range(len(args)) ]))
    return sum_pairwise

def wrap_sum_of_pairwise_with_coefficients(*args):
    """ Returns function that is sum of args functions. 

    Inputs:
        A list where the first half is a list of numbers and the second
    half is a list of functions that each depend on a distance (r).

    Returns:
        A function that evaluates the weighted sum of the inputted functions.
    The weights are the numbers in the first half of the inputted list.
    The returned function is also depends on the distance (r).
    NOT USED
    """
    def sum_pairwise_with_coefficients(r):
        return sum(np.array([ args[i]*args[i+(len(args)/2)](r) for i in range(len(args)/2) ]))
    return sum_pairwise_with_coefficients

def LJ12(r,r0):
    x = r0/r
    V = x**12
    return V

def LJ12_deriv(r,r0):
    x = r0/r
    V = (-12./r0)*(x**13)
    return V

def LJ1210(r,r0):
    x = r0/r
    V = 5.*(x**12) - 6.*(x**10)
    return V

def LJ1210_deriv(r,r0):
    x = r0/r
    V = (-60./r0)*((x**13) - (x**11))
    return V

def LJ1210rep(r,r0):
    x = r0/r
    V = np.zeros(x.shape,float)
    V[x > 1] = (5.*(x[x > 1]**12) - 6.*(x[x > 1]**10) + 2.)
    V[x <= 1] = -(5.*(x[x <= 1]**12) - 6.*(x[x <= 1]**10))
    return V

def LJ1210rep_deriv(r,r0):
    x = r0/r
    V = np.zeros(x.shape,float)
    V[x > 1] = -60.*(1./r0)*(x[x > 1]**13 - x[x > 1]**11)
    V[x <= 1] = 60.*(1./r0)*(x[x <= 1]**13 - x[x <= 1]**11)
    return V

def Gaussian(r,r0,width):
    V = -np.exp(-((r - r0)**2)/(2.*(width**2)))
    return V

def Gaussian_deriv(r,r0,width):
    V = ((r - r0)/(width**2))*np.exp(-((r - r0)**2)/(2.*(width**2)))
    return V
