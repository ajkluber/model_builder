''' Make pairwise_params strings for

Description:
    Construct pairwise_params file for compound interactions, such
as an LJ12 wall with Gaussian (Heiko's contacts). 

given structure (pdb file) and
set of contacts.
'''

import argparse
import pdb_parser

def get_LJ12_window_Gaussian_pairwise(pdb,contacts):
    ''' '''
    pairwise_distances = pdb_parser.get_pairwise_distances(pdb,contacts)
    model_param = 0
    ## Loop over contacts and create pair
    param_file_string = "#   i   j   param int_type  other_params\n"
    for i in range(len(contacts)):
        ## Each contact gets the interactions:
        ##  - LJ12_window_Gaussian      Int.Type = 8
        ##  - Gaussian                  Int.Type = 4
        param_file_string += "%5d%5d%5d%5d%s\n" % (i_idx,j_idx,model_param,int_type,other_param_string) 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--pdb', type=str, required=True, help='PDB file with structure.')
    parser.add_argument('--contacts', type=str, required=True, help='List of contacts.')
    args = parser.parse_args()

    pdbfile = args.pdb
    contactfile = args.contacts

    contacts = np.loadtxt(contactfile,dtype=int)
    if contacts.shape[0,:] == 4:
        contacts = contacts[:,1::2]

    pdb = pdb_parser.get_clean_CA(pdbfile)

    get_LJ12_window_Gaussian_pairwise(pdb,contacts)
