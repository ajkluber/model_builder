"""smog_AA_model

Description:
------------
    Reads smog.top files needed to run the smog all-atom protein model
in Gromacs. 

Example Usage:
--------------
    See project_tools/examples

"""

import numpy as np
import os
import shutil
from glob import glob

#import bead_representation as bead
#import bonded_potentials as bond
import smog_AA_pairwise_potentials as pairwise
import pdb_parser

#global SKIP_INTERACTIONS
#SKIP_INTERACTIONS = [1,8,9]

class smog_AA_model(object):
    """Interface for creating the gromacs files needed to run a coarse-grain
    simulation (.top,.gro,.ndx,table files). Also contains each term of the
    potential energy as a function for use in parameter fitting.
    """
    def __init__(self,**kwargs):
        self.path = os.getcwd()
        # Set any keyword argument given as an attribute. 
        for key in kwargs.iterkeys():
            setattr(self,key.lower(),kwargs[key])

        # Set remaining values to None
        need_to_define = ["disulfides","defaults","n_native_pairs","n_short_native_pairs","starting_gro",
                          "epsilon_bar","exclusions","contact_type",
                          "fitting_params_file","fitting_params",
                          "long_fitting_params","initial_t_array",
                          ]
        
        for thing in need_to_define:
            if not hasattr(self,thing):
                if thing == "exclusions":
                    self.exclusions = []
                else:
                    setattr(self,thing,None)
                    
        bonds_file = np.loadtxt(self.name+'/smog_files/smog_bonds.top')
        angles_file = np.loadtxt(self.name+'/smog_files/smog_angles.top')
        dihedrals_proper_file = np.loadtxt(self.name+'/smog_files/smog_dihedrals_proper.top')
        dihedrals_improper_file = np.loadtxt(self.name+'/smog_files/smog_dihedrals_improper.top')
        short_pairs_file = np.loadtxt(self.name+'/smog_files/smog_pairs_s.top') #short range pairs |i-j|<=8                  
#        fixed_long_pairs_file_1 = np.loadtxt(self.name+'/smog_files/smog_pairs_f1.top') #fixed part of long range pairs type 6
#        fixed_long_pairs_file_2 = np.loadtxt(self.name+'/smog_files/smog_pairs_f2.top') #removal of Gaussian part from above file (type 5)
        long_pairs_file = np.loadtxt(self.long_pairs_file) #variable (Gaussian) part of long-range pairs, |i-j|>8                                                                                  

############################                                                                                                                                     
# Here we make a decision not to optimize local contacts, and to just optimize long range ones.                                                                 
# One way of doing this is modifying the well width for Gaussian contacts (not being used)                                                                      
# Another way is to modify the value of the epsilons (the default way)                                                                                          
# Finally, long range contacts will be allowed to become repulsive (or the well diminished??)                                                                   
###########################                                                                                                                         
        self.using_sbm_gmx = True

        #Bonds parameters:                                                                                                                                      
        #func  r0(nm)  K_b                                                                                                                                                
        self.bonds = bonds_file[:,:2]
        self.bonds_params = bonds_file[:,2:]

        #Angles parameters:                                                                                                                                      
        #func  Theta_0 (deg)  K_a                                                                                                      
        self.angles = angles_file[:,:3]
        self.angles_params = angles_file[:,3:]

       #Dihedrals parameters:                                                                                                                                    
       #func   Phi_0   K_d   Mult.(propers)                                                                                             # ON HOLD FOR THE TIME BEING                                  
        self.dihedrals_proper = dihedrals_proper_file[:,:4]
        self.dihedrals_proper_params = dihedrals_proper_file[:,4:]

        self.dihedrals_improper = dihedrals_improper_file[:,:4]
        self.dihedrals_improper_params = dihedrals_improper_file[:,4:]
        # Pairs parameters:                                                                                                                                      
        # func                                                                                                                                                   
        # eps_c = intensity of contact (depth of well)                                                                                                           
        # r_0 = position of the minimum of the attractive well                                                                                                  
        # sigma = width of the attractive well                                                                                                                   
        # r_nc = position of the excluded volume hard wall                                                                
                                     
        self.short_pairs = short_pairs_file[:,:2]
        self.short_pairs_params = short_pairs_file[:,2:]
        self.long_pairs = long_pairs_file[:,:2]
        self.long_pairs_params = long_pairs_file[:,2:]

        self.exclusions = np.vstack((self.short_pairs,self.long_pairs))

        # Count pairs
        self.n_short_pairs = len(self.short_pairs)
        self.n_long_pairs = len(self.long_pairs)

        # Get pairwise_params (for long range contacts)
        self.n_long_fitting_params = len(self.long_pairs)

        n_res = int(open(self.name+'/smog_files/smog_atoms.top','r').readlines()[-1].split()[2])
        self.n_residues = n_res
    
        self.qref = np.zeros((self.n_residues,self.n_residues))
        
        residue_pairs_file = np.loadtxt('long_residue_contacts.dat')
        # Load model param interactions (
        for i in range(len(self.long_pairs)):
            a = int(residue_pairs_file[i][0])-1
            b = int(residue_pairs_file[i][1])-1
            
            if a > b:
                self.qref[a][b] = 1
            else:
                self.qref[b][a] = 1

        np.savetxt(self.name+"/smog_files/Qref_cryst.dat",self.qref,fmt="%4d",delimiter=" ")
#        np.savetxt(self.name+"/Qref_cryst.dat",self.qref,fmt="%4d",delimiter=" ")
        # Starting .gro file
        self.starting_gro = "smog.gro"

        #Generate fitting_params (pair index) for long_pairs
        self.long_fitting_params = range(self.n_long_pairs)
        self.long_pairwise_param_assignment = np.arange(self.n_long_pairs)
        
        # Generate list of pairwise_type (model_param_values)
        self.long_pairwise_type = self.long_pairs_params[:,0]
        # Epsilon
        self.long_model_param_values = self.long_pairs_params[:,1]
        # Sigma (equilibrium distance)
        self.long_pairwise_other_params = self.long_pairs_params[:,2]
        # Other Gaussian parameters
        self.long_pairwise_well_width = self.long_pairs_params[:,3]
#        self.long_pairwise_excluded_vol = fixed_long_pairs_file_1[:,-1]
        
        self.n_long_model_param = len(self.long_model_param_values)
      
        self.update_model_param_values(self.long_model_param_values)
        
        # Clean pdb for mutations
        self.pdb = self.name +'.pdb'
        self.cleanpdb_full = pdb_parser.get_clean_full(self.pdb)
        
    def update_model_param_values(self,new_model_param_values):
        """ If parameter changed sign, change the pairwise interaction type """
        # Switching between different interaction function types
        potential_type_switch = {2:3,3:2,5:9,9:5}
    
        # Loop over fitting_params only 
        for i in range(self.n_long_fitting_params):
            p_idx = self.long_fitting_params[i]
            if new_model_param_values[i] < 0.:
                    # If parameter changes sign then the pairwise_type is flipped.
                self.long_pairwise_type[p_idx] = potential_type_switch[self.long_pairwise_type[p_idx]]
            else:
                pass
                # Model parameters are always positive
            self.long_model_param_values[p_idx] = abs(new_model_param_values[i])   

        # Refresh everything that depends on model parameters
        self._determine_tabled_interactions()
        self._set_nonbonded_interactions()

    def save_simulation_files(self,savetables=True):
        """ Write all needed simulation files. """
        cwd = os.getcwd()
        relative_path = cwd.split("%s/" % self.path)[1]
        self.topology_file_location = "%s/smog_pairs_long" % relative_path
        smog_files = glob('{0}/{1}/smog_files/*'.format(self.path,self.name))                                              
        for item in smog_files:                                                                                                 
            shutil.copy(item,cwd)   
        os.remove("smog_pairs_l.top")
        os.remove("smog_pairs_long")
        os.remove("smog_bonds_rep.top")
        if os.path.isfile('frame.gro'):
            shutil.move('smog.gro','crystal.gro')
            shutil.move('frame.gro','smog.gro')
        open("smog_pairs_long","w").write(self.long_pairs_file_string)
        open("smog_pairs_l.top","w").write(self.long_pairs_top_string)
        open("smog_bonds_rep.top","w").write(self.long_bonds_rep_string)

        if savetables:
            # Save needed table files                                               
            np.savetxt("table.xvg",self.tablep,fmt="%16.15e",delimiter=" ")
            np.savetxt("tablep.xvg",self.tablep,fmt="%16.15e",delimiter=" ")
            for i in range(self.n_tables):
                np.savetxt(self.tablenames[i],self.tables[i],fmt="%16.15e",delimiter=" ")
                    

    def _determine_tabled_interactions(self):
        """Determine which interactions need to use a tableb_##.xvg file"""
        int_type = self.long_pairwise_type
        dont_table = [1,2,4,6,8,5,10]

        flag = np.zeros(len(int_type))
        for i in range(len(dont_table)):
            flag += (int_type == dont_table[i]).astype(int)

        self.tabled_interactions = np.zeros(self.n_long_pairs,float)
        for tbl_indx in (np.where(flag != 1))[0]:
            self.tabled_interactions[tbl_indx] = 1
        self.tabled_pairs = np.where(self.tabled_interactions == 1)[0]
        self.n_tables = len(self.tabled_pairs)

    def _set_nonbonded_interactions(self):
            """ Set all interaction functions """

        # File to save interaction parameters for pairwise potentials.                                                                                       
            pair_eps = []
            rep_pair_counter = 0

            self.long_pair_V = []
            self.long_pair_dV = []
            self.long_pairs_file_string = ""
            self.long_pairs_top_string = ""
            self.long_bonds_rep_string = ""

            for i in range(self.n_long_pairs):
                i_idx = int(self.long_pairs[i][0])
                j_idx = int(self.long_pairs[i][1])
                model_param = self.long_pairwise_param_assignment[i]
                int_type = int(self.long_pairwise_type[i])
                eps = self.long_model_param_values[i]
                sigma = self.long_pairwise_other_params[i]
                well_width = self.long_pairwise_well_width[i]
#                excluded_vol = self.long_pairwise_excluded_vol[i]
                    
                self.long_pairs_file_string += '{0:4d}   {1:4d}    {2:2d}    {3:2.12e}   {4:2.12e}   {5:2.12e}\n'.format(i_idx,j_idx,int_type,eps,sigma,well_width)
            # Wrap the pairwise potentials so that only distance needs to be input.                         
                if int_type == 9:
                    rep_pair_counter+=1
                    self.long_bonds_rep_string = '     {0:>4d}     {1:>4d}     {2:1d}            {3:4d} {4:2.9e}\n'.format(i_idx,j_idx,int_type,rep_pair_counter,eps)
                else:
                    self.long_pairs_top_string += '{0:4d}   {1:4d}    {2:2d}    {3:2.12e}   {4:2.12e}   {5:2.12e}\n'.format(i_idx,j_idx,int_type,eps,sigma,well_width)

                self.long_pair_V.append(pairwise.wrap_pairwise(pairwise.get_pair_potential(int_type),\
                                                               sigma,well_width))
                self.long_pair_dV.append(pairwise.wrap_pairwise(pairwise.get_pair_potential_deriv(int_type),\
                                            sigma,well_width))

            # Assign pairwise interaction strength from model parameters                                                                                     
                pair_eps.append(self.long_model_param_values[self.long_pairwise_param_assignment[i]])
                
            self.long_pair_eps = np.array(pair_eps)

 # File to save model parameters                                                                                                                      
##### TO-DO : GENERATE TANG INTERACTION TABLES           
            self._generate_interaction_tables()
            
    def _generate_interaction_tables(self):
        """ Generates tables of user-defined potentials                                                                                                                                                                              
        """
        self.tablep = self._get_LJ1210_table()
        self.LJtable = self.tablep
        r = np.arange(0,20.0,0.002)
        self.tables = []
        self.tablenames = []
        for i in range(self.n_tables):
            pair_indx = self.tabled_pairs[i]
            table_name = "table_b%d.xvg" % (i+1)
            self.tablenames.append(table_name)

            table = np.zeros((len(r),3),float)
            table[1:,0] = r[1:]
            table[10:,1] = self.long_pair_V[pair_indx](r[10:])
            table[10:,2] = -1.*self.long_pair_dV[pair_indx](r[10:])
            self.tables.append(table)

    def _get_LJ1210_table(self):
        """ LJ1210 interaction table """
        r = np.arange(0.0,100.0,0.002)
        r[0] = 1
        table = np.zeros((len(r),7),float)
        table[:,0] = r
        table[:,1] = 1./r
        table[:,2] = 1./(r**2)
        table[:,3] = -1./(r**10)
        table[:,4] = -10./(r**11)
        table[:,5] = 1./(r**12)
        table[:,6] = 12./(r**13)
        table[:5,1:] = 0
        table[0,0] = 0
        return table

    def _get_tabled_string(self):

        """ Generate the topology files to specify table interactions. """
        # Add special nonbonded table interactions.                                                                             
 
        tabled_string = "; tabled interactions pairs below\n"
        for i in range(self.n_tables):
            pair_indx = self.tabled_pairs[i]
            pair = self.long_pairs[pair_indx]
            eps = self.long_pair_eps[pair_indx]
            tabled_string += "{0:4d}   {1:4d}    {2:2d}    {3:2.12e}   {4:2.12e}\n".format(pair[0],pair[1],9,i+1,eps)
        return tabled_string
     
    def save_table_files(self):

        """Write only table files"""
        np.savetxt("table.xvg",self.tablep,fmt="%16.15e",delimiter=" ")
        np.savetxt("tablep.xvg",self.tablep,fmt="%16.15e",delimiter=" ")
        for i in range(self.n_tables):
            np.savetxt(self.tablenames[i],self.tables[i],fmt="%16.15e",delimiter=" ")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build a SMOG model.')
    parser.add_argument('--name', type=str, required=True, help='pdb')
    parser.add_argument('--pairs', type=str, required=True, help='pairs')
    args = parser.parse_args() 

    name = args.name
#    pdb = "%s.pdb" % name
#    pairsfile = args.pairs



    model = smog_AA_model()
    #model.save_simulation_files()
