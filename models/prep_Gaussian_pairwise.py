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

#
#    count = 0
#    C = np.zeros((len(C_nn),len(C_nn)))*np.nan
#    for i in range(len(C_nn)):
#        atm_i = i + 1
#        for j in range(i+1,len(C_nn)):
#            atm_j = j + 1
#            flag = sum((pairs[:,0] == atm_i).astype(int)*(pairs[:,1] == atm_j).astype(int))
#            native_dist = np.linalg.norm(atm_coords[i] - atm_coords[j])
#            if (flag == 0) and (native_dist > nn_cutoff):
#                C[j,i] = eps_random[count]
#                count += 1

def add_nonnative_interactions(args,pdb,native_pairs,pairwise_param_file_string,model_param_file_string,model_param,random=False):

    nn_cutoff = 0.85

    pdb_info = pdb_parser.get_coords_atoms_residues(pdb)
    atm_coords = pdb_info[0]
    res_types = pdb_info[4]
    eps_nn = []
    # First determine the distribution of non-native contact energies.
    for i in range(len(res_types)):
        atm_i = i + 1
        for j in range(i+1,len(res_types)):
            atm_j = j + 1
            flag = sum((native_pairs[:,0] == atm_i).astype(int)*(native_pairs[:,1] == atm_j).astype(int))
            native_dist = np.linalg.norm(atm_coords[i] - atm_coords[j])
            if (flag == 0) and (native_dist > nn_cutoff):
                res_i = res_types[i]
                res_j = res_types[j]
                if random:
                    gamma = np.random.normal(loc=args.nonnative_mean,scale=np.sqrt(args.nonnative_var))
                else:
                    gamma = rp.get_awsem_direct_contact_gamma(res_i,res_j)
                eps_nn.append(gamma)

    # Shift and scale the non-native contact energies to control the mean and
    # variance.
    eps_nn = np.array(eps_nn)
    eps_nn_avg = np.mean(eps_nn)
    eps_nn_var = np.var(eps_nn)
    eps_nn -= eps_nn_avg
    eps_nn *= np.sqrt(args.nonnative_var/eps_nn_var)
    eps_nn += args.nonnative_mean

    rNC_nn = 0.4
    r0_nn = rNC_nn + 0.1
    width0 = 0.075
    count = 0
    C_nn = np.zeros((len(res_types),len(res_types)))*np.nan
    # Create lines for pairwise_params to specify non-native interactions.
    for i in range(len(res_types)):
        atm_i = i + 1
        for j in range(i+1,len(res_types)):
            atm_j = j + 1
            flag = sum((native_pairs[:,0] == atm_i).astype(int)*(native_pairs[:,1] == atm_j).astype(int))
            native_dist = np.linalg.norm(atm_coords[i] - atm_coords[j])
            if (flag == 0) and (native_dist > nn_cutoff):
                gamma = eps_nn[count]
                C_nn[j,i] = gamma
                # compound_LJ12_Gaussian takes rNC, r0, width0
                LJ12_Gaussian_other_params = "%10.5f%10.5f%10.5f" % (rNC_nn,r0_nn,width0)
                pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (atm_i,atm_j,model_param,8,LJ12_Gaussian_other_params) 
                model_param += 1
                model_param_file_string += "%10.5f\n" % 1.

                if gamma < 0.:
                    interaction_type = 5
                else:
                    interaction_type = 4
                eps = abs(gamma)
                Gaussian_other_params = "%10.5f%10.5f" % (r0_nn,width0)
                pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (atm_i,atm_j,model_param,interaction_type,Gaussian_other_params) 
                model_param_file_string += "%10.5f\n" % eps
                model_param += 1
                count += 1

    return pairwise_param_file_string, model_param_file_string, eps_nn, C_nn


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

def get_compound_LJ12_Gaussian_pairwise(pdb,pairs,model_param=0):
    """ """

    width0 = 0.05
    pairwise_distances = pdb_parser.get_pairwise_distances(pdb,pairs)
    model_param_value = 1.
    pairs_rNC = 0.4*np.ones(len(pairs))

    # Loop over pairs and create pair
    pairwise_param_file_string = "#   i   j   param int_type  other_params\n"
    model_param_file_string = "# model parameters\n"
    for i in range(len(pairs)):
        # Each contact gets the interactions:
        #  - compound_LJ12_Gaussian      Int.Type = 8
        #  - Gaussian                    Int.Type = 4
        i_idx = pairs[i][0]
        j_idx = pairs[i][1]
        r0 = pairwise_distances[i] 
        rNC = pairs_rNC[i]

        # compound_LJ12_Gaussian takes rNC, r0, width0
        LJ12_Gaussian_other_params = "%10.5f%10.5f%10.5f" % (rNC,r0,width0)
        pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (i_idx,j_idx,model_param,8,LJ12_Gaussian_other_params) 
        model_param += 1
        model_param_file_string += "%10.5f\n" % model_param_value

        Gaussian_other_params = "%10.5f%10.5f" % (r0,width0)
        pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (i_idx,j_idx,model_param,4,Gaussian_other_params) 
        model_param_file_string += "%10.5f\n" % model_param_value
        model_param += 1

    return pairwise_param_file_string, model_param_file_string, model_param

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--pdb', type=str, required=True, help='PDB file with structure.')
    parser.add_argument('--pairs', type=str, required=True, help='List of pairs.')
    parser.add_argument('--noclean',action='store_true',default=False,help='Do not clean the pdb.')
    parser.add_argument('--nonnative', action='store_true', default=False, help='Optional. Add flavored non-native interactions.')
    parser.add_argument('--nonnative_mean', type=float, default=None, help='Optional. Non-native interaction mean.')
    parser.add_argument('--nonnative_var',  type=float, default=None, help='Optional. Non-native interaction std dev.')
    args = parser.parse_args()

    pdbfile = args.pdb
    pairsfile = args.pairs
    name = pdbfile.split(".pdb")[0]

    pairs = np.loadtxt(pairsfile,dtype=int)
    if pairs.shape[1] == 4:
        pairs = pairs[:,1::2]

    if not args.noclean:
        pdb = pdb_parser.get_clean_CA(pdbfile)
    else:
        pdb = open(pdbfile,"r").read()

    pairwise_param_file_string, model_param_file_string, model_param = get_compound_LJ12_Gaussian_pairwise(pdb,pairs)

    if args.nonnative:
        if (args.nonnative_mean is None) or (args.nonnative_var is None):
            raise IOError("Need to specify --nonnative_mean --nonnative_var for nonnative interactions.")
        pairwise_param_file_string, model_param_file_string, eps_nn, C_nn = add_nonnative_interactions(args,pdb,pairs,pairwise_param_file_string,model_param_file_string,model_param,random=True)

    if args.nonnative:
        open("%s_pairwise_params_%.2f_%.2f" % (name,args.nonnative_mean,args.nonnative_var),"w").write(pairwise_param_file_string)
        open("%s_model_params_%.2f_%.2f" % (name,args.nonnative_mean,args.nonnative_var),"w").write(model_param_file_string)
    else:
        open("%s_pairwise_params" % name,"w").write(pairwise_param_file_string)
        open("%s_model_params" % name,"w").write(model_param_file_string)

