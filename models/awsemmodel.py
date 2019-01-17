from __future__ import print_function
import numpy as np
import os
import shutil
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
        self.param_source = param_path
        
    def add_fragment_memory(self, param_path, mem_file, max_frag_length=9, cycle = True, fragment_memory_scale=0.1):
        """ Add the fragment memory terms to the model from a .mem file
        
        Parameters
        ----------
        param_path : str
            Path to directory containing the .mem file
        mem_file : str
            Name of the .mem file.
        
        """
        print("Adding memory file %s from %s" % (mem_file, param_path))
        self.Hamiltonian.fragment_memory_scale = fragment_memory_scale
        cwd = os.getcwd()
        os.chdir(param_path)
        f = open(mem_file, "r")
        f.readline() # read the [memories] line
        self.fragment_info = []
        for line in f: #cycle through until you're out of lines
            info = line.strip().split() 
            self.fragment_info.append(info)
            mem_pdb = info[0]
            ##ASSUMPTION: First index is protein, second is for fragment
            #subtract one to convert to python indices
            protein_index = int(info[1]) - 1 
            frag_index = int(info[2]) - 1
            length = int(info[3])
            weight = float(info[4])
            traj = md.load(mem_pdb)
            if length > max_frag_length and cycle:
                print("Cycling Through the sequence, exceeding max_frag_length")
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
                self.Hamiltonian.add_fragment_memory(traj, protein_index, frag_index, length, weight)
        os.chdir(cwd)
    
    def copy_awsem_parameters(self, new_directory, start_directory=None):
        if start_directory is None:
            start_directory = self.param_source
            
        if not os.path.isdir(new_directory):
            os.mkdir(new_directory)
            
        param_files = os.listdir(start_directory)
        for pfile in param_files:
            pfile_str = "%s/%s" % (start_directory, pfile)
            if os.path.isfile(pfile_str):
                shutil.copy(pfile_str, new_directory)
                
        self.param_source = new_directory
    
    def write_new_fragment_memory(self, new_param_path, name):
        frag_str = "[Memories]\n"
        for idx, info in enumerate(self.fragment_info):
            for jdx in range(4):
                frag_str += info[jdx]
                frag_str += " "
            frag_str += "%.6f\n" % self.Hamiltonian.fragment_potentials[idx].weight
        f = open("%s.mem" % name, "w")
        f.write(frag_str)
        
    def write_new_gammas(self, new_param_path):
        gamma_str = self._get_gammas_file_format()
        f = open("%s/gamma.dat"%new_param_path, "w")
        f.write(gamma_str)
        f.close()
                
    def _get_gammas_file_format(self):
        gamma_direct = self.Hamiltonian.gamma_direct
        ##debug
        ##gamma_direct=np.zeros((20,20))
        ##
        gamma_str = ""
        for i in range(20):
            for j in range(i,20):
                gval = gamma_direct[i,j]
                gamma_str += "%8.5f  %8.5f\n" % (gval, gval)
        
        gamma_str += "\n"
        gamma_water = self.Hamiltonian.gamma_water
        gamma_protein = self.Hamiltonian.gamma_protein
        
        for i in range(20):
            for j in range(i,20):
                wval = gamma_water[i,j]
                pval = gamma_protein[i,j]
                gamma_str += "%8.5f  %8.5f\n" % (pval, wval)
        
        return gamma_str
