
import numpy as np
import models

def _get_pairwise_params():
    pair_params = np.loadtxt("pairwise_params",skiprows=1)
    model_params = np.loadtxt("model_params",skiprows=1)
    
    contacts = pair_params[:,:2].astype(int)
    pairwise_param_assignment = pair_params[:,2].astype(int)
    pairwise_type = pair_params[:,3].astype(int)
    pairwise_other_params = pair_params[:,4:]
    
    return contacts,pairwise_param_assignment,model_params,pairwise_type,pairwise_other_params

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build a structure-based model.')
    parser.add_argument('--name', type=str, required=True, help='pdb')
    args = parser.parse_args() 

    name = args.name
    pdb = "%s.pdb" % name

    contacts,pairwise_param_assignment,model_params,pairwise_type,pairwise_other_params = _get_pairwise_params()

    model = models.SmogCalpha.SmogCalpha(pdb=pdb,pairwise_other_params=pairwise_other_params,contacts=contacts,model_param_values=model_params)
    model.save_simulation_files()
