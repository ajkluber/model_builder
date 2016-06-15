import argparse
import numpy as np

import mdtraj as md

from model_builder.models.mappings import CalphaMapping

def native_gaussians(pairs, rNC, r0, width, model_param):
    """Create pairwise params string describe native Gaussian contacts

    Parameters
    ----------
    pairs : np.ndarray int (n_pairs, 2)
        Indices of pairs of coarse-grain atoms that are in native contact.
    rNC : float
        Excluded volume radius for native contacts.
    r0 : np.ndarray float (n_pairs,)
        Native contact distances.
    width : float
        Width of of the Gaussian well for the native contact. 
    model_param : int
        Starting index of the model parameter (epsilon).

    """
    # Add native contacts
    eps_native = 1 
    eps_string = "# epsilon\n" 
    params_string = "#   i   j   param int_type  other_params\n" 
    for i in range(len(pairs)):
        params_string += "{:>5d}{:>5d}{:>7d}{:>15s}{:>10.5f}{:>10.5f}{:>10.5f}\n".format(
                                    pairs[i,0], pairs[i,1], model_param, 
                                    "LJ12GAUSSIAN", rNC, r0[i], width)
        eps_string += "{:>10.5f}\n".format(eps_native)
        model_param += 1
    return params_string, eps_string, model_param

def add_nonnative_interactions(ref_traj, n_residues, pairs, std_dev_eps_nn,
                    rNC_nn, r0_nn, width_nn, model_param=0, params_string="",
                    eps_string="", nn_cutoff=0.85, eps_NN=0):
    """Create pairwise params string describe non-native interactions

    The strength of non-native contacts is drawn from a normal distribution

    Parameters
    ----------
    ref_traj : obj. mdtraj.Trajectory
        Reference structure 
    n_residues : int
        Number of residues
    pairs : np.ndarray int (n_pairs, 2)
        Indices of pairs of coarse-grain atoms that are in native contact.
    std_dev_eps_nn : float
        The standard deviation of the non-native contact distribution. 
    rNC_nn : float
        Excluded volume radius for non-native contacts.
    r0_nn : np.ndarray float (n_pairs,)
        Non-native contact distances.
    width_nn : float
        Width of of the Gaussian well for the non-native contact. 
    nn_cutoff : float
        Non-native contacts are excluded if they are within a cutoff distance
        in the native state. This is that cutoff.
    model_param : int, optional
        Starting index of the model parameter (epsilon).
    params_string : str, optional
        The pairwise parameters string.
    model_string : str, optional
        The model parameters string.
    eps_NN : float, optional
        The center of the epsilon distibution.

    """

    # tuples work better for comparisons
    tup_pairs = [ tuple(pair) for pair in pairs ]

    nonnative_pairs = []
    for i in range(n_residues):
        for j in range(i + 4, n_residues):
            if not (tuple([i + 1, j + 1]) in tup_pairs):
                # Only non-native pair interactions
                r0_nn_ref = np.linalg.norm(ref_traj.xyz[0, i, :] - ref_traj.xyz[0, j, :])

                if (r0_nn_ref <= nn_cutoff):
                    # exclude non-native contacts that are too close in the native
                    # state.
                    continue
                else:
                    # Non-native strength is a random variable
                    eps_nn = np.random.normal(loc=eps_NN, scale=std_dev_eps_nn)
                    nonnative_pairs.append([i + 1, j + 1, eps_nn])
                    if eps_nn < 0:
                        # Repulsive tanh interaction 
                        params_string += "{:>5d}{:>5d}{:>7d}{:>15s}{:>10.5f}{:>10.5f}{:>10.5f}\n".format(
                                                i + 1, j + 1, model_param, "LJ12TANHREP", 
                                                rNC_nn, r0_nn, width_nn)
                    else:
                        # Attractive Gaussian interaction 
                        params_string += "{:>5d}{:>5d}{:>7d}{:>15s}{:>10.5f}{:>10.5f}{:>10.5f}\n".format(
                                                i + 1, j + 1, model_param, "LJ12GAUSSIAN",
                                                rNC_nn, r0_nn, width_nn)
                    eps_string += "{:>10.5f}\n".format(abs(eps_nn))
                    model_param += 1

    return params_string, eps_string, model_param, nonnative_pairs

def make_pairwise_params(path_to_pdb, path_to_contacts, std_dev_eps_nn):
    """Create pairwise params file with non-native contacts


    Parameters
    ----------
    path_to_pdb : str
        The path to the pdb file 
    path_to_contacts : str
        The path to the contacts file 
    std_dev_eps_nn : float
        The standard deviation of the non-native contact distribution. 

    """

    # calculate native contact distances
    ref_pdb = md.load(path_to_pdb)
    mapping = CalphaMapping(ref_pdb.top)
    ref_traj = mapping.map_traj(ref_pdb)
    
    pairs = np.loadtxt(path_to_contacts, dtype=int)
    if pairs.shape[1] == 4:
        pairs = pairs[:,(1,3)] 

    r0 = md.compute_distances(ref_traj, pairs - 1)[0]

    model_param = 0

    # Add Native interactions as Gaussians
    rNC = 0.4       # excluded volume radius
    width = 0.05    # width of the Gaussian well
    params_string, eps_string, model_param = native_gaussians(pairs, rNC, r0, width, model_param)

    # Add non-native interactions to the pairwise_params file
    rNC_nn = 0.4            # excluded volume radius
    r0_nn = rNC_nn + 0.1    # non-native contact distance
    width_nn = 0.075        # width of gaussian/tanh well
    nn_cutoff = 0.85        # cutoff radius to exclude near-native non-native contacts

    nonnative_pairs = []
    if std_dev_eps_nn != 0: 
        params_string, eps_string, model_param, nonnative_pairs = \
            add_nonnative_interactions(ref_traj, mapping.top.n_residues, pairs,
                        std_dev_eps_nn, rNC_nn, r0_nn, width_nn, model_param=model_param,
                        params_string=params_string, eps_string=eps_string,
                        nn_cutoff=nn_cutoff) 

    return params_string, eps_string, nonnative_pairs

if __name__ == "__main__":
    """Format the pairwise_params file for Gaussian contacts"""
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb")
    parser.add_argument("contacts")
    parser.add_argument("b_value")
    args = parser.parse_args()
    
    params_string, eps_string, nonnative_pairs = make_pairwise_params(args.pdb, args.contacts, float(args.b_value))

    name = args.pdb.split(".pdb")[0]
    with open("{}_pairwise_params".format(name), "w") as fout:
        fout.write(params_string)

    with open("{}_model_params".format(name), "w") as fout:
        fout.write(eps_string)
