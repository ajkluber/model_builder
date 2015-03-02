""" Make pairwise_params strings for

Description:
    Construct pairwise_params file for compound interactions, such
as an LJ12 wall with Gaussian (Heiko's pairs). 

Inputs:
    pdb file string
    pairs

Outputs:
    pairwise_params file string

TODO:
    Allow for specifying the model_param of each interaction.

"""

import argparse
import numpy as np

import pdb_parser
import residue_properties as rp

def get_Gaussian_pairwise(pdb,pairs,model_param=0):
    """ """
    width0 = 0.05
    pairwise_distances = pdb_parser.get_pairwise_distances(pdb,pairs)
    model_param_value = 1.

    ## Loop over pairs and create pair
    pairwise_param_file_string = "#   i   j   param int_type  other_params\n"
    model_param_file_string = "# model parameters\n"
    for i in range(len(pairs)):
        ## Each contact gets the interactions:
        ##  - Gaussian                    Int.Type = 4
        i_idx = pairs[i][0]
        j_idx = pairs[i][1]
        r0 = pairwise_distances[i] 

        Gaussian_other_params = "%10.5f%10.5f" % (r0,width0)
        pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (i_idx,j_idx,model_param,4,Gaussian_other_params) 
        model_param_file_string += "%10.5f\n" % model_param_value
        model_param += 1

    return pairwise_param_file_string, model_param_file_string

def get_compound_LJ12_Gaussian_pairwise(pdb,pairs,model_param=0,native_sterics=False):
    """ """
    width0 = 0.05
    pairwise_distances = pdb_parser.get_pairwise_distances(pdb,pairs)
    model_param_value = 1.
    if native_sterics:
        res_types = pdb_parser.get_coords_atoms_residues(pdb)[4]
        pairs_rNC = np.zeros(len(pairs))
        for i in range(len(pairs)):
            i_idx = pairs[i][0]
            j_idx = pairs[i][1]
            radii_i = rp.residue_CB_radii(res_types[i_idx - 1])
            radii_j = rp.residue_CB_radii(res_types[j_idx - 1])
            pairs_rNC = radii_i + radii_j
    else:
        pairs_rNC = 0.4*np.ones(len(pairs))

    ## Loop over pairs and create pair
    pairwise_param_file_string = "#   i   j   param int_type  other_params\n"
    model_param_file_string = "# model parameters\n"
    for i in range(len(pairs)):
        ## Each contact gets the interactions:
        ##  - compound_LJ12_Gaussian      Int.Type = 8
        ##  - Gaussian                    Int.Type = 4
        i_idx = pairs[i][0]
        j_idx = pairs[i][1]
        r0 = pairwise_distances[i] 
        rNC = pairs_NC[i]

        ## compound_LJ12_Gaussian takes rNC, r0, width0
        LJ12_Gaussian_other_params = "%10.5f%10.5f%10.5f" % (rNC,r0,width0)
        pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (i_idx,j_idx,model_param,8,LJ12_Gaussian_other_params) 
        model_param += 1
        model_param_file_string += "%10.5f\n" % model_param_value

        Gaussian_other_params = "%10.5f%10.5f" % (r0,width0)
        pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (i_idx,j_idx,model_param,4,Gaussian_other_params) 
        model_param_file_string += "%10.5f\n" % model_param_value
        model_param += 1


    return pairwise_param_file_string, model_param_file_string

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--pdb', type=str, required=True, help='PDB file with structure.')
    parser.add_argument('--pairs', type=str, required=True, help='List of pairs.')
    parser.add_argument('--bead_repr', type=str, required=True, help='Bead representation.')
    args = parser.parse_args()

    pdbfile = args.pdb
    pairsfile = args.pairs
    bead_repr = args.bead_repr
    name = pdbfile.split(".pdb")[0]

    pairs = np.loadtxt(pairsfile,dtype=int)
    if pairs.shape[1] == 4:
        pairs = pairs[:,1::2]

    if bead_repr == "CA":
        pdb = pdb_parser.get_clean_CA(pdbfile)
    elif bead_repr == "CACB":
        pdb = pdb_parser.get_clean_CA_center_of_mass_CB(pdbfile)
    else:
        raise IOErorr("--bead_repr must be CA or CACB")

    pairwise_param_file_string, model_param_file_string = get_compound_LJ12_Gaussian_pairwise(pdb,pairs)

    open("%s_pairwise_params" % name,"w").write(pairwise_param_file_string)
    open("%s_model_params" % name,"w").write(model_param_file_string)

