"""Utilities for bonded potential terms"""

import numpy as np

def distance(coords,i_idx,j_idx):
    dist = np.linalg.norm(coords[i_idx] - coords[j_idx])
    return dist

def angle(coords,i_idx,j_idx,k_idx):
    """ Compute the angle between 3 atoms in degrees"""
    xkj = coords[k_idx] - coords[j_idx]
    xkj /= np.linalg.norm(xkj)
    xij = coords[i_idx] - coords[j_idx]
    xij /= np.linalg.norm(xij)
    theta = (180./np.pi)*np.arccos(np.dot(xkj, xij))
    return theta

def dihedral(coords,i_idx,j_idx,k_idx,l_idx):
    """ Compute the dihedral angle between planes formed by 4 atoms in degrees """
    v21 = coords[j_idx] - coords[i_idx]
    v31 = coords[k_idx] - coords[i_idx]
    v32 = coords[k_idx] - coords[j_idx]
    v42 = coords[l_idx] - coords[j_idx]
    v21xv31 = np.cross(v21,v31)
    v21xv31 /= np.linalg.norm(v21xv31)
    v32xv42 = np.cross(v32,v42)
    v32xv42 /= np.linalg.norm(v32xv42)
    if np.dot(v21,v32xv42) < 0.:
        sign = -1.
    else:
        sign = 1.
    phi = 180. + sign*(180./np.pi)*np.arccos(np.dot(v21xv31,v32xv42))
    return phi
