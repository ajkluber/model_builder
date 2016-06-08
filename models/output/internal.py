import os
import numpy as np

import mdtraj as md

class InternalFiles(object):
    
    def __init__(self, model):
        self.model = model
    
    def write_pairwise_parameters(self, suffix=""):
        pairwise_string = self._generate_pairwise_string()
        model_params_string = self._generate_model_params_string()
        
        f = open("pairwise_params%s" % suffix, "w")
        f.write(pairwise_string)
        f.close()
        
        f = open("model_params%s" % suffix, "w")
        f.write(model_params_string)
        f.close()
        
        
    def _generate_pairwise_string(self):
        pairwise_string = "#    pairs         param         potential_type   other_params\n"
        count = 0
        for pot in self.model.Hamiltonian._pairs:
            index_i = pot.atmi.index + 1
            index_j = pot.atmj.index + 1
            potential_type = pot.prefix_label
            other_parameters = pot.other_params
            
            current_str = "%6d  %6d  %12d %18s     " % (index_i, index_j, count, potential_type)
            
            for param in other_parameters:
                current_str += "%.6f  " % param
            
            current_str += "\n"
            
            pairwise_string += current_str
            count += 1
        
        return pairwise_string
    
    def _generate_model_params_string(self):
        model_params_string = "# model parameters\n"
        for eps in self.model.Hamiltonian._epsilons:
            model_params_string += "%f\n" % eps
            
        return model_params_string
