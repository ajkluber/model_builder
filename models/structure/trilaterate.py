import numpy as np


def trilaterate(v1, v2, v3, d1, d2, d3, left_handed=True):
    """Find position of fourth atom with three distances of three other atoms.

    Parameters
    ----------
    v1 : np.ndarray(3) 
        Coordinates of first atom. Used as origin.
    v2 : np.ndarray(3) 
        Coordinates of second atom. Aligned to x-axis.
    v3 : np.ndarray(3) 
        Coordinates of third atom. Placed in the xy-plane.
    d1 : float
        Distance to first atom.
    d2 : float
        Distance to second atom.
    d3 : float
        Distance to third atom.

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Trilateration
    """

    # define orthonormal basis vectors e_x, e_y, e_z
    v2_x = np.linalg.norm(v2 - v1)
    e_x = (v2 - v1)/v2_x

    v3_x = np.dot(e_x, v3 - v1)

    e_y = (v3 - v1 - v3_x*e_x)/np.linalg.norm(v3 - v1 - v3_x*e_x)

    v3_y = np.dot(e_y, v3 - v1)

    e_z = np.cross(e_x, e_y)

    # Should be unit norm.
    #print np.linalg.norm(e_x), np.linalg.norm(e_y), np.linalg.norm(e_z) 

    v4_x = (d1**2 - d2**2 + v2_x**2)/(2.*v2_x)
    v4_y = ((d1**2 - d3**2 + v3_x**2 + v3_y**2)/(2.*v3_y)) - (v3_x/v3_y)*v4_x
    v4_z = np.sqrt(d1**2 - v4_x**2 - v4_y**2)

    # two roots indicate left- and right-handed chiralities. 
    if left_handed: 
        # amino acid sidechains are left-handed.
        v4 = v1 + v4_x*e_x + v4_y*e_y - v4_z*e_z
    else:
        v4 = v1 + v4_x*e_x + v4_y*e_y + v4_z*e_z

    return v4

if __name__ == "__main__":
    import mdtraj as md
    import matplotlib.pyplot as plt

    traj = md.load("complex.pdb")

    # Testing that this algorithm reproduces the CB positions of all
    # non-glycine residues in a test structure.
    all_d1 = []
    all_d2 = []
    all_d3 = []
    for res in traj.top.residues:

        if res.name != "GLY":
            C_idx = [ atom.index for atom in res.atoms if (atom.name == "C") ][0]
            N_idx = [ atom.index for atom in res.atoms if (atom.name == "N") ][0]
            CA_idx = [ atom.index for atom in res.atoms if (atom.name == "CA") ][0]
            CB_idx = [ atom.index for atom in res.atoms if (atom.name == "CB") ][0]

            v1 = traj.xyz[0,N_idx,:]
            v2 = traj.xyz[0,CA_idx,:]
            v3 = traj.xyz[0,C_idx,:]
            v4_true = traj.xyz[0,CB_idx,:]

            # distance constraints to backbone atoms
            d1 = np.linalg.norm(v4_true - v1)
            d2 = np.linalg.norm(v4_true - v2)
            d3 = np.linalg.norm(v4_true - v3)

            # v4 is position of a fourth atom.
            v4 = trilaterate(v1, v2, v3, d1, d2, d3)

            #print np.allclose(v4, v4_true), np.linalg.norm(v4 - v4_true)
            
            all_d1.append(d1)
            all_d2.append(d2)
            all_d3.append(d3)

    
    # The average values
    d1 = 0.2470955
    d2 = 0.1533931
    d3 = 0.2510052
