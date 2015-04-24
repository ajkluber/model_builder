import numpy as np

import pdb_parser

def get_smog_pairs(pdb,pairs,pairwise_strengths,distances=None,rNC=0.04,sigma=0.05):
    """Get [pairs] section for SBM gmx hack

    Parameters
    ----------
    pdb : str
        String in PDB file format that has contents of one protein chain. 
        Assumed to be one chain, numbering starting at 1. Can also be a 
        filename that ends in .pdb.
    pairs : array 
        Array of shape (n_pairs,2) that contains the atom indices to calculate
        the distance between.
    rNC : array or float
        Array of shape (n_pairs) that contains the non-contact radius for each 
        interaction. If a float is passed it will be used .
    rNC : array or float
        Array of shape (n_pairs) that contains the non-contact radius for each 
        interaction. If a float is passed it will be used .

    """
    if isinstance(rNC,float):
        rNC = rNC*np.ones(pairs.shape[0])
    else:
        if len(rNC) != pairs.shape[0]:
            raise IOError("rNC argument must be array with the same size as the number of pairs.")

    if isinstance(sigma,float):
        sigma = sigma*np.ones(pairs.shape[0])
    else:
        if len(sigma) != pairs.shape[0]:
            raise IOError("sigma argument must be array with the same size as the number of pairs.")

    if distances is None:
        distances = pdb_parser.get_pairwise_distances(pdb,pairs)
    else:
        if len(distances) != pairs.shape[0]:
            raise IOError("distances argument must be array with the same size as the number of pairs.")

    ftype = 6
    smog_string = "[ pairs ]\n"
    for i in xrange(len(pairs)):
        smog_string += "%5d%5d%5d %8.5f %8.5f %8.5f %8.5e\n" % \
                        (pairs[i][0],pairs[i][1],ftype,pairwise_strengths[i],distances[i],sigma[i],rNC[i]**12)

    return smog_string

if __name__ == "__main__":
    name = "SH3"
    pdb = open("%s.pdb" % name,"r").read()
    pairs = np.loadtxt("%scaca_cbcb_pairs" % name, dtype=int)
    strengths = np.ones(len(pairs))

    smog_string = get_smog_pairs(pdb,pairs,strengths)

    with open("smog_pairs","w") as fout:
        fout.write(smog_string)
