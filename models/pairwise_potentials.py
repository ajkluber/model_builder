""" Utilities for nonbonded/pairwise potential terms.

Potentials to add: 
- Repulsive tanh function by Ryan Cheng to complement Gaussian  DONE
- Attractive and repulsive LJ126 potentials DONE
- Compound functions i.e. Heiko's product functions  DONE
- Compound (product) functions with Gaussian barrier. 

- Add second derivatives

"""

import numpy as np

##############################################################################
# Utilities and wrappers
##############################################################################

def get_pair_potential(code):
    """ Returns pairwise potential function"""
    potential = {1:LJ12,
                 2:LJ1210,
                 3:LJ1210rep,
                 4:Gaussian,
                 5:Cheng_rep,
                 6:LJ126,
                 7:LJ126rep,
                 8:compound_LJ12_Gaussian,
                 9:compound_LJ12_double_Gaussian,
                 10:compound_double_Gaussian,
                 11:FRET_Efficiency}[code]
    return potential

def get_pair_potential_deriv(code):
    """ Returns derivative of pairwise potential function"""
    potential = {1:LJ12_deriv,
                 2:LJ1210_deriv,
                 3:LJ1210rep_deriv,
                 4:Gaussian_deriv,
                 5:Cheng_rep_deriv,
                 6:LJ126_deriv,
                 7:LJ126rep_deriv,
                 8:compound_LJ12_Gaussian_deriv,
                 9:compound_LJ12_double_Gaussian_deriv,
                 10:compound_double_Gaussian_deriv,
                 11:FRET_Efficiency_deriv}[code]
    return potential

def get_switched_pair_potential(code):
    """ Returns the "switched" pairwise potential function. NOT DONE/USED"""
    potential = {1:LJ12,2:LJ1210,3:LJ1210rep,4:Gaussian,5:Cheng_rep}[code]
    return potential

def wrap_pairwise(old_pairwise,*args):
    """Wrap pairwise function so only r needs to be passed 

    Parameters:
    -----------
    old_pairwise : function
        A function that takes one distance (r) and a variable number of
        positional arguments depends on (such as equilibrium distance or width
        of well).
    args : tuple
        A tuple that contains values of position arguments for old_pairwise

    Returns:
    --------
    new_pairwise : function
        A function that only depends on the distance, with the remaining 
        positional parameters hard coded.
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

##############################################################################
# Pairwise interaction potentials
##############################################################################

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

def Cheng_rep(r,r0,width):
    alpha = 1./width      
    r0prime = r0 + width
    V = 0.5*(np.tanh(-alpha*(r - r0prime)) + 1.)
    return V

def Cheng_rep_deriv(r,r0,width):      
    alpha = 1./width      
    r0prime = r0 + width
    V = -0.5*alpha*(1. - (np.tanh(-alpha*(r - r0prime)))**2)
    return V

def LJ126(r,r0):
    x = r0/r
    V = (x**12) - 2.*(x**6)
    return V

def LJ126_deriv(r,r0):
    x = r0/r
    V = (-12./r0)*((x**13) - (x**7))
    return V

def LJ126rep(r,r0):
    """ ONLY PLACEHOLDER. NOT DONE """
    x = r0/r
    V = (x**12) - 2.*(x**6)
    return V

def LJ126rep_deriv(r,r0):
    """ ONLY PLACEHOLDER. NOT DONE """
    x = r0/r
    V = (-12./r0)*((x**13) - (x**7))
    return V


##############################################################################
# Compound pairwise potentials. 
##############################################################################

def compound_LJ12_Gaussian(r,rNC,r0,width0):
    return LJ12(r,rNC)*(1. + Gaussian(r,r0,width0))

def compound_LJ12_Gaussian_deriv(r,rNC,r0,width0):
    return LJ12_deriv(r,rNC)*(1. + Gaussian(r,r0,width0)) + LJ12(r,rNC)*Gaussian_deriv(r,r0,width0)

def compound_LJ12_double_Gaussian(r,rNC,r0,width0,r1,width1):
    # Not tested in simulation
    return LJ12(r,rNC)*(1. + Gaussian(r,r0,width0))*(1. + Gaussian(r,r1,width1))

def compound_LJ12_double_Gaussian_deriv(r,rNC,r0,width0,r1,width1):
    return LJ12_deriv(r,rNC)*(1. + Gaussian(r,r0,width0))*(1. + Gaussian(r,r1,width1)) + \
           LJ12(r,rNC)*Gaussian_deriv(r,r0,width0)*(1. + Gaussian(r,r1,width1)) + \
           LJ12(r,rNC)*(1. + Gaussian(r,r0,width0))*Gaussian_deriv(r,r1,width1)
v
def compound_double_Gaussian(r,rNC,r0,width0,r1,width1):
    return Gaussian(r,r0,width0)*(1. + Gaussian(r,r1,width1))

def compound_double_Gaussian_deriv(r,rNC,r0,width0,r1,width1):
    return Gaussian_deriv(r,r0,width0)*(1. + Gaussian(r,r1,width1)) +\
           Gaussian(r,r0,width0)*Gaussian_deriv(r,r1,width1)

##############################################################################
# Special pairwise potentials. 
##############################################################################

def FRET_Efficiency(r,r0):
    x = r/r0
    V = 1./(1. + x**6)
    return V

def FRET_Efficiency_deriv(r,r0):
    x = r/r0
    V = (-6./r0)*(x**5)/((1. + x**6)**2)
    return V

def FRET_Efficiency_moment(r,r0,eff_mean,moment):
    V = (FRET_Efficiency(r,r0) - eff_mean)**moment
    return V

def FRET_Efficiency_moment_deriv(r,r0,eff_mean,moment):
    V = moment*((FRET_Efficiency(r,r0) - eff_mean)**(moment - 1))*FRET_Efficiency_deriv(r,r0)
    return V
