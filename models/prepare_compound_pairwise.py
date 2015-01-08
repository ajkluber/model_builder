''' Make pairwise_params strings for

Description:
    Construct pairwise_params file for compound interactions, such
as an LJ12 wall with Gaussian (Heiko's contacts). 

Inputs:
    pdb file string
    contacts

Outputs:
    pairwise_params file string

TODO:
    Allow for specifying the model_param of each interaction.

'''

import argparse
import numpy as np

import pdb_parser

def get_compound_LJ12_Gaussian_pairwise(pdb,contacts,model_param=0):
    ''' '''
    rNC = 0.4
    width0 = 0.05
    pairwise_distances = pdb_parser.get_pairwise_distances(pdb,contacts)
    model_param_value = 1.

    ## Loop over contacts and create pair
    pairwise_param_file_string = "#   i   j   param int_type  other_params\n"
    model_param_file_string = "# model parameters\n"
    for i in range(len(contacts)):
        ## Each contact gets the interactions:
        ##  - compound_LJ12_Gaussian      Int.Type = 8
        ##  - Gaussian                    Int.Type = 4
        i_idx = contacts[i][0]
        j_idx = contacts[i][1]
        r0 = pairwise_distances[i] 

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
    parser.add_argument('--contacts', type=str, required=True, help='List of contacts.')
    args = parser.parse_args()

    pdbfile = args.pdb
    contactfile = args.contacts
    name = pdbfile.split(".pdb")[0]

    contacts = np.loadtxt(contactfile,dtype=int)
    if contacts.shape[0] == 4:
        contacts = contacts[:,1::2]

    pdb = pdb_parser.get_clean_CA(pdbfile)

    pairwise_param_file_string, model_param_file_string = get_compound_LJ12_Gaussian_pairwise(pdb,contacts)

    open("%s_pairwise_params" % name,"w").write(pairwise_param_file_string)
    open("%s_model_params" % name,"w").write(model_param_file_string)

