import numpy as np

import mdtraj as md

from model_builder.models.mappings.calpha import CalphaMapping

if __name__ == "__main__":
    int_type = "LJ12GAUSSIAN"
    rNC = 0.4
    width = 0.05
    
    # Calculate the native contact distances
    ref_pdb = md.load("SH3.pdb")
    mapping = CalphaMapping(ref_pdb.top)
    
    pairs = np.loadtxt("SH3.contacts", dtype=int)
    n_pairs = len(pairs)

    r0 = md.compute_distances(mapping.map_traj(ref_pdb), pairs - 1)[0]

    # Format the pairwise_params file for Gaussian contacts
    params_string = "#   i   j   param int_type  other_params\n" 
    for i in range(n_pairs):
        params_string += "{:>5d}{:>5d}{:>7d}{:>15s}{:>10.5f}{:>10.5f}{:>10.5f}\n".format(
                                    pairs[i,0], pairs[i,1], i, int_type, rNC, r0[i], width)

    with open("SH3_pairwise_params", "w") as fout:
        fout.write(params_string)
