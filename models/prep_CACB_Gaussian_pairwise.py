import argparse
import numpy as np

import residue_properties as rp
import pdb_parser

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of directory.')
    parser.add_argument('--n_caca_pairs', type=int, required=True, help='Num .')
    parser.add_argument('--cushion', type=float, default=0.2, help='Name of directory.')
    args = parser.parse_args()

    name = args.name
    n_caca_pairs = args.n_caca_pairs
    cushion = args.cushion
    avg_radius = rp.residue_radii["AVERAGE"]

    pdb = "%s.pdb" % name
    pairs = np.loadtxt("%scaca_cbcb_pairs" % name,dtype=int)

    cacbpdb = pdb_parser.get_clean_CA_center_of_mass_CB(pdb)
    pairwise_distances = pdb_parser.get_pairwise_distances(cacbpdb,pairs)

    width0 = 0.05
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
        if (i + 1) <= n_caca_pairs:
            rNC = 0.2*2
        else:
            rNC = 2.*avg_radius

        if r0 < (1. + cushion)*rNC:
            # Scale down excluded volume of close contacts
            rNC /= (1. + cushion)

        # compound_LJ12_Gaussian takes rNC, r0, width0
        LJ12_Gaussian_other_params = "%10.5f%10.5f%10.5f" % (rNC,r0,width0)
        pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (i_idx,j_idx,model_param,8,LJ12_Gaussian_other_params) 
        model_param += 1
        model_param_file_string += "%10.5f\n" % model_param_value

        Gaussian_other_params = "%10.5f%10.5f" % (r0,width0)
        pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (i_idx,j_idx,model_param,4,Gaussian_other_params) 
        model_param_file_string += "%10.5f\n" % model_param_value
        model_param += 1

    open("%s_pairwise_params" % name,"w").write(pairwise_param_file_string)
    open("%s_model_params" % name,"w").write(model_param_file_string)
