import argparse
import numpy as np

import residue_properties as rp
import pdb_parser

def add_nonnative_interactions(args,pairwise_param_file_string,model_param_file_string,model_param,pdb,cacbpdb):
    # Add non-native flavored interactions
    avg_radius = rp.residue_radii["AVERAGE"]
    rNC_nonnative = 2.*avg_radius
    r0_nonnative = rNC_nonnative + 0.1
    width0 = 0.075

    # Non-native interactions between CB-CB pairs.
    pdb_info = pdb_parser.get_coords_atoms_residues(cacbpdb)
    atm_coords = pdb_info[0]
    atm_indxs = pdb_info[1]
    atm_types = pdb_info[2]
    res_indxs = pdb_info[3]
    res_types = pdb_info[4]
    CA_indxs = atm_indxs[atm_types == "CA"]
    CB_indxs = atm_indxs[atm_types == "CB"]
    n_residues = len(np.unique(np.array(res_indxs,copy=True)))
    n_atoms = len(atm_indxs)
    n_CB = len(CB_indxs)
    n_CA = len(CA_indxs)
    C = np.zeros((n_residues,n_residues),float)  
    nonnative_gamma = []
    for i in range(n_CB):
        for j in range(i+1,n_CB):
            atm_idx_i = CB_indxs[i]
            atm_idx_j = CB_indxs[j]
            idx_i = np.where(atm_indxs == atm_idx_i)[0][0]
            idx_j = np.where(atm_indxs == atm_idx_j)[0][0]

            # Non-native atom pair
            res_i = res_types[idx_i]
            res_j = res_types[idx_j]

            if np.linalg.norm(atm_coords[idx_i - 1] - atm_coords[idx_j - 1]) > 0.8: 
                # Exclude non-native interactions that may alter native state energy.
                gamma = rp.get_awsem_direct_contact_gamma(res_i,res_j)
                C[res_indxs[idx_j] - 1,res_indxs[idx_i] - 1] = gamma
                nonnative_gamma.append(gamma)

                # Decide on sum scaling for non-native contacts.

                if args.nonnative:
                    # compound_LJ12_Gaussian takes rNC, r0, width0
                    LJ12_Gaussian_other_params = "%10.5f%10.5f%10.5f" % (rNC_nonnative,r0_nonnative,width0)
                    pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (atm_idx_i,atm_idx_j,model_param,8,LJ12_Gaussian_other_params) 
                    model_param += 1
                    model_param_file_string += "%10.5f\n" % 1.

                    if gamma < 0.:
                        interaction_type = 5
                    else:
                        interaction_type = 4
                    eps = abs(gamma)
                    Gaussian_other_params = "%10.5f%10.5f" % (r0_nonnative,width0)
                    pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (atm_idx_i,atm_idx_j,model_param,interaction_type,Gaussian_other_params) 
                    model_param_file_string += "%10.5f\n" % eps
                    model_param += 1

    # Add contacts from glycine CA to all other CB's.
    rNC_gly_nn = 0.4
    r0_gly_nn = rNC_gly_nn + 0.1
    gly_atm_idxs = atm_indxs[res_types == "GLY"]
    gly_res_idxs = res_indxs[res_types == "GLY"]
    CB_res_types = res_types[atm_types == "CB"]
    CB_res_idxs = res_indxs[atm_types == "CB"]
    for i in range(len(gly_atm_idxs)):
        idx_i = gly_atm_idxs[i]
        res_idx_i = gly_res_idxs[i]
        for j in range(n_CB):
            res_idx_j = CB_res_idxs[j]
            res_j = CB_res_types[j]
            idx_j = CB_indxs[j]

            gamma = rp.get_awsem_direct_contact_gamma("GLY",res_j)
            if idx_i > idx_j: 
                atm_idx_i = idx_j
                atm_idx_j = idx_i
                temp = res_idx_i
                res_idx_i = res_idx_j
                res_idx_j = temp
            else:
                atm_idx_i = idx_i
                atm_idx_j = idx_j

            if np.linalg.norm(atm_coords[idx_i - 1] - atm_coords[idx_j - 1]) > 0.8: 
                if args.nonnative:
                    # compound_LJ12_Gaussian takes rNC, r0, width0
                    LJ12_Gaussian_other_params = "%10.5f%10.5f%10.5f" % (rNC_gly_nn,r0_gly_nn,width0)
                    pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (atm_idx_i,atm_idx_j,model_param,8,LJ12_Gaussian_other_params) 
                    model_param += 1
                    model_param_file_string += "%10.5f\n" % 1.
                    C[res_idx_j - 1,res_idx_i - 1] = gamma

                    if gamma < 0.:
                        interaction_type = 5
                    else:
                        interaction_type = 4
                    eps = abs(gamma)
                    Gaussian_other_params = "%10.5f%10.5f" % (r0_gly_nn,width0)
                    pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (atm_idx_i,atm_idx_j,model_param,interaction_type,Gaussian_other_params) 
                    model_param_file_string += "%10.5f\n" % eps
                    model_param += 1

    # Add interactions between glycines
    for i in range(len(gly_atm_idxs)):
        idx_i = gly_atm_idxs[i]
        atm_idx_i = idx_i
        res_idx_i = gly_res_idxs[i]
        for j in range(i+1,len(gly_atm_idxs)):
            idx_j = gly_atm_idxs[j]
            atm_idx_i = idx_j
            res_idx_j = gly_res_idxs[i]
            gamma = rp.get_awsem_direct_contact_gamma("GLY","GLY")
            if np.linalg.norm(atm_coords[idx_i - 1] - atm_coords[idx_j - 1]) > 0.8: 
                if args.nonnative:
                    # compound_LJ12_Gaussian takes rNC, r0, width0
                    LJ12_Gaussian_other_params = "%10.5f%10.5f%10.5f" % (rNC_gly_nn,r0_gly_nn,width0)
                    pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (atm_idx_i,atm_idx_j,model_param,8,LJ12_Gaussian_other_params) 
                    model_param += 1
                    model_param_file_string += "%10.5f\n" % 1.
                    C[res_idx_j - 1,res_idx_i - 1] = gamma

                    if gamma < 0.:
                        interaction_type = 5
                    else:
                        interaction_type = 4
                    eps = abs(gamma)
                    Gaussian_other_params = "%10.5f%10.5f" % (r0_gly_nn,width0)
                    pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (atm_idx_i,atm_idx_j,model_param,interaction_type,Gaussian_other_params) 
                    model_param_file_string += "%10.5f\n" % eps
                    model_param += 1

    return pairwise_param_file_string, model_param_file_string, C

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of pdb.')
    parser.add_argument('--pairsfile', type=str, default=None, help='Optional. Name of pairs file.')
    parser.add_argument('--nonnative', action='store_true', default=False, help='Optional. Add flavored non-native interactions.')
    parser.add_argument('--nonnative_scaling', type=float,default=0.01, help='Optional. Scaling constant for non-native interactions.')
    #parser.add_argument('--cushion', type=float, default=0.2, help='Name of directory.')
    args = parser.parse_args()

    name = args.name

    #cushion = args.cushion

    pdb = "%s.pdb" % name

    # Get CACB pdb.
    cacbpdb = pdb_parser.get_clean_CA_center_of_mass_CB(pdb)

    if args.pairsfile is not None:
        # User defined CACA CBCB contact file
        #pairs = np.loadtxt("%scaca_cbcb_pairs" % name,dtype=int)
        pairs = np.loadtxt("%s" % name,dtype=int)
    else:
        # Determine contacts by cutoff map
        pairs = pdb_parser.get_CACB_contacts_cutoff(pdb)

    # Contact contact distances
    pairwise_distances = pdb_parser.get_pairwise_distances(cacbpdb,pairs)

    # Create Gaussian interactions
    width0 = 0.05
    rNC = 0.4
    model_param_value = 1.
    model_param = 0
    pairwise_param_file_string = "#   i   j   param int_type  other_params\n"
    model_param_file_string = "# model parameters\n"
    for i in range(len(pairs)):
        # Each contact gets the interactions:
        #  - compound_LJ12_Gaussian      Int.Type = 8
        #  - Gaussian                    Int.Type = 4
        i_idx = pairs[i,0]
        j_idx = pairs[i,1]
        r0 = pairwise_distances[i] 

        # compound_LJ12_Gaussian takes rNC, r0, width0
        LJ12_Gaussian_other_params = "%10.5f%10.5f%10.5f" % (rNC,r0,width0)
        pairwise_param_file_string += "%5d%5d%7d%5d%s\n" % (i_idx,j_idx,model_param,8,LJ12_Gaussian_other_params) 
        model_param += 1
        model_param_file_string += "%10.5f\n" % 1.

        Gaussian_other_params = "%10.5f%10.5f" % (r0,width0)
        pairwise_param_file_string += "%5d%5d%7d%5d%s\n" % (i_idx,j_idx,model_param,4,Gaussian_other_params) 
        model_param_file_string += "%10.5f\n" % model_param_value
        model_param += 1

    pairwise_param_file_string, model_param_file_string, C = add_nonnative_interactions(args,pairwise_param_file_string,model_param_file_string,model_param,pdb,cacbpdb)

    if args.nonnative:
        open("%s_%.2f_pairwise_params" % (name,args.nonnative_scaling),"w").write(pairwise_param_file_string)
        open("%s_%.2f_model_params" % (name,args.nonnative_scaling),"w").write(model_param_file_string)
    else:
        open("%s_pairwise_params" % name,"w").write(pairwise_param_file_string)
        open("%s_model_params" % name,"w").write(model_param_file_string)
