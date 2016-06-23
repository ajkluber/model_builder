import numpy as np
import os
import mdtraj as md

import potentials 

from model import Model            


class AwsemModel(Model):

    def __init__(self, topology):
        Model.__init__(self, topology, bead_repr="AWSEM")
        self.Hamiltonian = potentials.AwsemHamiltonian(self.mapping.top)
        self.Hamiltonian._set_charged_residues(self.mapping._charged_residues)

    def source_parameters(self, param_path, parameterize=True): 
        # read in parameters from source directory
        self.Hamiltonian._source_parameters(param_path)

    def add_fragment_memory(self, param_path, mem_file, max_frag_length=9, cycle = True, fragment_memory_scale=0.1):
        """ Add the fragment memory terms to the model from a .mem file
        
        Parameters
        ----------
        param_path : str
            Path to directory containing the .mem file
        mem_file : str
            Name of the .mem file.
        
        """
        
        self.Hamiltonian.fragment_memory_scale = fragment_memory_scale
        cwd = os.getcwd()
        os.chdir(param_path)
        f = open(mem_file, "r")
        f.readline() # read the [memories] line
        for line in f: #cycle through until you're out of lines
            info = line.strip().split() 
            mem_pdb = info[0]
            ##ASSUMPTION: First index is protein, second is for fragment
            #subtract one to convert to python indices
            protein_index = int(info[1]) - 1 
            frag_index = int(info[2]) - 1
            length = int(info[3])
            weight = float(info[4])
            traj = md.load(mem_pdb)
            if length > max_frag_length and cycle:
                print "Sliding Through Both Sequences"
                total_length = length
                length = max_frag_length
                go = True
                count = 0
                while go:
                    p_idx = protein_index + count
                    f_idx = frag_index + count
                    self.Hamiltonian.add_fragment_memory(traj, p_idx, f_idx, length, weight)
                    if count+max_frag_length == total_length:
                        go = False
                    count += 3
            else:   
                if length > max_frag_length:
                    print "No Cycling, exceeding max_frag_length" 
                self.Hamiltonian.add_fragment_memory(traj, protein_index, frag_index, length, weight)
        os.chdir(cwd)
         
        
