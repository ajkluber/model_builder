
import numpy as np
import models

def get_pairwise_params(pairwise_params_file,model_params_file):
    pair_params = np.loadtxt(pairwise_params_file)
    model_param_values = np.loadtxt(model_params_file)
    
    pairs = pair_params[:,:2].astype(int)
    pairwise_param_assignment = pair_params[:,2].astype(int)
    pairwise_type = pair_params[:,3].astype(int)
    pairwise_other_parameters = pair_params[:,4:]
    
    return pairs,pairwise_param_assignment,model_param_values,pairwise_type,pairwise_other_parameters

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

    model = models.SmogCalpha.SmogCalpha(pdb=pdb,pairs=pairs,
                    pairwise_param_assignment=pairwise_param_assignment,
                    pairwise_type=pairwise_type,
                    pairwise_other_parameters=pairwise_other_parameters,
                    model_param_values=model_param_values)

    model.save_simulation_files()
