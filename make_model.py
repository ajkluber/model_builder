
import numpy as np
import models

def get_pairwise_params(pairwise_params_file,model_params_file):
    """ Grab pairwise_params from file. """
    model_param_values = np.loadtxt(model_params_file)

    p_lines = [ x.rstrip("\n") for x in open(pairwise_params_file,"r").readlines() ]

    pairs = []
    pairwise_param_assignment = []
    pairwise_type = [] 
    pairwise_other_params = []

    for p in p_lines[1:]:
        data = p.split() 
        pairs.append([int(data[0]),int(data[1])])
        pairwise_param_assignment.append(int(data[2]))
        pairwise_type.append(int(data[3]))
        temp = []
        for otherparam in data[4:]:
            temp.append(float(otherparam))
        pairwise_other_params.append(temp)

    pairs = np.array(pairs) 
    pairwise_param_assignment = np.array(pairwise_param_assignment)
    pairwise_type = np.array(pairwise_type)

    return pairs,pairwise_param_assignment,model_param_values,pairwise_type,pairwise_other_params

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build a structure-based model.')
    parser.add_argument('--name', type=str, required=True, help='pdb')
    args = parser.parse_args() 

    name = args.name
    pdb = "%s.pdb" % name
    pairwise_params_file = "pairwise_params"
    model_params_file = "model_params"

    pairs,pairwise_param_assignment,model_param_values,pairwise_type,pairwise_other_parameters = get_pairwise_params(pairwise_params_file,model_params_file)

    model = models.CoarseGrainedModel.CoarseGrainedModel(pdb=pdb,pairs=pairs,
                    pairwise_param_assignment=pairwise_param_assignment,
                    pairwise_type=pairwise_type,
                    pairwise_other_parameters=pairwise_other_parameters,
                    model_param_values=model_param_values,bead_repr="CACB")

    model.save_simulation_files()
