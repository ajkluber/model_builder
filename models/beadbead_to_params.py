""" 
November 2014
For backwards compatability, split all BeadBead files into 

interaction_params
model_params

"""
import os
import glob

import numpy as np


def convert_beadbead(filename):

    LJtype_to_int_type = {1:2,-1:3}
    beadbead = np.loadtxt(filename,dtype=str)
    params = beadbead[:,4]
    use = params != "ss"
    params = params[use].astype(int)
    pairs = beadbead[use,:2].astype(int)
    epsilons = beadbead[use,6].astype(float)
    sigmas = beadbead[use,5].astype(float)

    LJtype = (beadbead[use,7].astype(float)).astype(int)

    int_params = "#   i   j   param int_type  other_params\n"
    model_params = "# model parameters\n"
    for i in range(sum(use.astype(int))):
        i_idx = pairs[i][0]
        j_idx = pairs[i][1]
        model_p = params[i]
        int_type = LJtype_to_int_type[LJtype[i]]
        other_param_string = " %10.5f" % sigmas[i]

        int_params += "%5d%5d%5d%5d%s\n" % (i_idx,j_idx,model_p,int_type,other_param_string)

        model_params += "%10.5f\n" % epsilons[i]

    return int_params, model_params

if __name__ == "__main__":
    ## Traverse directory directory tree and convert all BeadBead.dat's into
    ## interation_params and model_params
    ## Find all BeadBead.dat files within three layers of current directory 
    ## and add the new file format alongside them.

    beads = glob.glob("*/*/*/BeadBead.dat") + glob.glob("*/*/BeadBead.dat") \
          + glob.glob("*/*/*/NewBeadBead*.dat") + glob.glob("*/*/NewBeadBead*.dat")

    #beads = ["1FMK/Tf_0/119_0/BeadBead.dat"]
    for n in range(len(beads)):
        fullpath = beads[n]
        dir = os.path.dirname(beads[n])
        int_params, model_params = convert_beadbead(fullpath)
        open(dir+"/pairwise_params","w").write(int_params)
        open(dir+"/model_params","w").write(model_params)

