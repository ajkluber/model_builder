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
import matplotlib.pyplot as plt

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

def add_random_nonnative_interactions(args,pdb,native_pairs,pairwise_param_file_string,model_param_file_string,model_param,random=False):

    nn_cutoff = 0.85

    pdb_info = pdb_parser.get_coords_atoms_residues(pdb)
    atm_coords = pdb_info[0]
    res_types = pdb_info[4]
    if not (args.var == 0.00):
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
                        gamma = np.random.normal(loc=args.mean,scale=np.sqrt(args.var))
                    else:
                        gamma = rp.get_awsem_direct_contact_gamma(res_i,res_j)
                    eps_nn.append(gamma)

        # Shift and scale the non-native contact energies to control the mean and
        # variance.
        eps_nn = np.array(eps_nn)
        eps_nn_avg = np.mean(eps_nn)
        eps_nn_var = np.var(eps_nn)
        eps_nn -= eps_nn_avg
        eps_nn *= np.sqrt(args.var/eps_nn_var)
        eps_nn += args.mean
    else:
        eps_nn = [0]*(len(res_types)**2)

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
                if args.var == 0.00:
                    gamma = 0
                    interaction_type = 4
                else:
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

def get_compound_LJ12_Gaussian_pairwise(args,pdb,pairs,model_param=0,cacb=False,flavored=False ):
    """ """

    width0 = 0.05
    pairwise_distances = pdb_parser.get_pairwise_distances(pdb,pairs)
    pdb_info = pdb_parser.get_coords_atoms_residues(pdb)
    n_residues = len(pdb_info[4]) 
    C_native = np.zeros((n_residues,n_residues),float)
    stuff = pdb_parser.get_coords_atoms_residues(pdb)
    atomtype = stuff[2]
    residuetype = stuff[4]
    if args.random_native:
        epsilons = np.random.normal(loc=args.avg_native_eps,scale=np.sqrt(args.var),size=len(pairs))
        epsilons *= (args.avg_native_eps/np.mean(epsilons))
        epsilons -= args.avg_native_eps
        epsilons *= np.sqrt(args.var)/np.std(epsilons)
        epsilons += args.avg_native_eps
    else:
        epsilons = np.ones(len(pairs),float)
        
    pairs_rNC = 0.4*np.ones(len(pairs))
    if cacb == True: #if cacb, set excluded volumne for flavored interactions or average interactions according to selection rule
        for i in range(len(pairs)):
            ith_idx = pairs[i][0] - 1
            jth_idx = pairs[i][1] - 1
            #set params for calpha
            ith_param = 0.28
            jth_param = 0.28 
            #if cb, change
            if atomtype[ith_idx] == "CB":
                if flavored == True:
                    ith_param = rp.residue_cacb_effective_interaction[residuetype[ith_idx]]
                else:
                    ith_param = rp.residue_cacb_effective_interaction["AVERAGE"]
            if atomtype[jth_idx] == "CB":
                if flavored == True:
                    jth_param = rp.residue_cacb_effective_interaction[residuetype[jth_idx]]
                else:
                    jth_param = rp.residue_cacb_effective_interaction["AVERAGE"]

            pairs_rNC[i] = (ith_param * jth_param) ** 0.5

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
        gamma = epsilons[i]
        
        C_native[j_idx - 1, i_idx - 1] = gamma

        # compound_LJ12_Gaussian takes rNC, r0, width0
        LJ12_Gaussian_other_params = "%10.5f%10.5f%10.5f" % (rNC,r0,width0)
        pairwise_param_file_string += "%5d%5d%7d%5d%s\n" % (i_idx,j_idx,model_param,8,LJ12_Gaussian_other_params) 
        model_param += 1
        model_param_file_string += "%10.5f\n" % 1.

        if gamma < 0.:
            interaction_type = 5
        else: 
            interaction_type = 4
        model_param_value = abs(gamma)

        Gaussian_other_params = "%10.5f%10.5f" % (r0,width0)
        pairwise_param_file_string += "%5d%5d%7d%5d%s\n" % (i_idx,j_idx,model_param,interaction_type,Gaussian_other_params) 
        model_param_file_string += "%10.5f\n" % model_param_value
        model_param += 1

    return pairwise_param_file_string, model_param_file_string, model_param, C_native

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--pdb', 
                        type=str, 
                        required=True, 
                        help='PDB file with structure.')

    parser.add_argument('--pairs', 
                        type=str, 
                        required=True, 
                        help='List of pairs.')
    
    parser.add_argument('--noclean',
                        action='store_true',
                        default=False,
                        help='Do not clean the pdb.')

    parser.add_argument('--random_native', 
                        action='store_true', 
                        default=False, 
                        help='Optional. Add random native heterogeneity.')

    parser.add_argument('--random_nonnative', 
                        action='store_true', 
                        default=False, 
                        help='Optional. Add random non-native interactions.')

    parser.add_argument('--avg_native_eps', 
                        type=float, 
                        default=1., 
                        help='Optional. Average native interaction strength.')
    
    parser.add_argument('--mean', 
                        type=float, 
                        default=None, 
                        help='Optional. Non-native interaction mean.')

    parser.add_argument('--var',  
                        type=float, 
                        default=None, 
                        help='Optional. Non-native interaction variance.')
    
    parser.add_argument('--cacb',
                        action='store_true',
                        default=False,
                        help='Making a CACB model')
                        
    parser.add_argument('--flavored',
                        action='store_true',
                        default=False,
                        help='Use flavored radii for CA-CB model')
    
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

    pairwise_param_file_string, model_param_file_string, model_param, C_native = get_compound_LJ12_Gaussian_pairwise(args,pdb,pairs, cacb=args.cacb, flavored=args.flavored)

    cmap = plt.get_cmap("jet")
    cmap.set_bad(color="gray",alpha=1.)
    if args.random_nonnative:
        if (args.mean is None) or (args.var is None):
            raise IOError("Need to specify --mean --var for random_nonnative interactions.")
        pairwise_param_file_string, model_param_file_string, eps_nn, C_nn = add_random_nonnative_interactions(args,pdb,pairs,pairwise_param_file_string,model_param_file_string,model_param,random=True)

        C_nn = np.ma.masked_invalid(C_nn)
        plt.figure()
        plt.pcolormesh(C_nn,cmap=cmap)
        plt.xlim(0,C_nn.shape[0])
        plt.ylim(0,C_nn.shape[0])
        plt.xlabel("Residue i")
        plt.ylabel("Residue j")
        cbar = plt.colorbar()
        cbar.set_label("Interaction strength")
        plt.title("%s random non-native  $\\overline{\\epsilon_{nn}} = %.2f$ $\\sigma_{\\epsilon_{nn}}^2 = %.2f$ " \
                    % (name,args.mean,args.var))
        plt.savefig("%s_random_%.2f_%.2f.png" % (name,args.mean,args.var),format="png")
        plt.savefig("%s_random_%.2f_%.2f.pdf" % (name,args.mean,args.var),format="pdf")
        plt.savefig("%s_random_%.2f_%.2f.eps" % (name,args.mean,args.var),format="eps")
        #plt.show()

    if args.random_native:
        plt.figure()
        plt.pcolormesh(C_native,cmap=cmap)
        plt.xlim(0,C_native.shape[0])
        plt.ylim(0,C_native.shape[0])
        plt.xlabel("Residue i")
        plt.ylabel("Residue j")
        cbar = plt.colorbar()
        cbar.set_label("Interaction strength")
        plt.title("%s random native  $\\overline{\\epsilon_{nat}} = %.2f$ $\\sigma_{\\epsilon_{nat}}^2 = %.2f$ " \
                    % (name,args.avg_native_eps,args.var))
        plt.savefig("%s_random_native_%.2f_%.2f.png" % (name,args.avg_native_eps,args.var),format="png")
        plt.savefig("%s_random_native_%.2f_%.2f.pdf" % (name,args.avg_native_eps,args.var),format="pdf")
        plt.savefig("%s_random_native_%.2f_%.2f.eps" % (name,args.avg_native_eps,args.var),format="eps")
        #plt.show()


    if args.random_nonnative:
        open("%s_pairwise_params_%.2f_%.2f" % (name,args.mean,args.var),"w").write(pairwise_param_file_string)
        open("%s_model_params_%.2f_%.2f" % (name,args.mean,args.var),"w").write(model_param_file_string)
    elif args.random_native:
        open("%s_pairwise_params_%.2f_%.2f" % (name,args.avg_native_eps,args.var),"w").write(pairwise_param_file_string)
        open("%s_model_params_%.2f_%.2f" % (name,args.avg_native_eps,args.var),"w").write(model_param_file_string)
    else:
        open("%s_pairwise_params" % name,"w").write(pairwise_param_file_string)
        open("%s_model_params" % name,"w").write(model_param_file_string)

