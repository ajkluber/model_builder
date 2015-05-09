"""CoarseGrainedModel

Description:
------------
    Generates topol.top needed to run coarse-grained protein model
in Gromacs.

Example Usage:
--------------
    See project_tools/examples

"""

import numpy as np
import os

import bead_representation as bead
import bonded_potentials as bond
import pairwise_potentials as pairwise
import pdb_parser

global SKIP_INTERACTIONS
SKIP_INTERACTIONS = [1,8,9]

class CoarseGrainedModel(object):
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
        need_to_define = ["disulfides","defaults","n_native_pairs","starting_gro",
                          "epsilon_bar","exclusions","contact_type",
                          "backbone_params","initial_T_array",
                          "fitting_params_file","fitting_params",
                          "pairwise_params_file_location","model_params_file_location"]

        for thing in need_to_define:
            if not hasattr(self,thing):
                if thing == "exclusions":
                    self.exclusions = []
                elif thing == "backbone_params":
                    # SMOG values Kb = 2E4; Ka = 40, Kd = 1
                    self.backbone_params = ["Kb","Ka","Kd"]
                    self.backbone_param_vals = {"Kb":20000.,"Ka":40.,"Kd":1} 
                else:
                    setattr(self,thing,None)

        # Set backbone parameters
        bead.set_bonded_interactions(self)

        # Generate topology file and table files.
        self.n_pairs = len(self.pairs)
        self._check_pair_opts()
        self._determine_native_pairs()
        self._determine_tabled_interactions()
        self._set_nonbonded_interactions()

    def _check_pair_opts(self):
        """ Set default pairwise interaction terms """

        # If not specified, allow all parameters to be fit.
        if self.fitting_params is None:
            self.fitting_params = range(self.n_pairs)
        self.n_fitting_params = len(self.fitting_params)

        # If defaults is False then all these need to be defined.
        if self.defaults:
            if self.verbose:
                print "  Using defaults: LJ1210 interactions, homogeneous strength 1, each param is free."
            self.pairwise_distances = pdb_parser.get_pairwise_distances(self.cleanpdb,self.pairs)
            self.model_param_values = np.ones(self.n_pairs,float)
            self.pairwise_type = 2*np.ones(self.n_pairs,int)
            self.pairwise_param_assignment = np.arange(self.n_pairs)     
            self.pairwise_other_parameters = [ tuple({self.pairwise_distances[x]}) for x in range(self.n_pairs) ]
        else:
            needed = [hasattr(self,"model_param_values"),hasattr(self,"pairwise_type"), \
                    hasattr(self,"pairwise_param_assignment"),hasattr(self,"pairwise_other_parameters")]
            if not all(needed):
                err  = "If not using defualts then the following need to be inputted:\n"
                err += "     model_param_values\n"
                err += "     pairwise_type\n"
                err += "     pairwise_param_assignment\n"
                err += "     pairwise_other_parameters"
                raise IOError(err)

        self.give_smaller_excluded_volume = []
        #if self.exclusions == []:
        for i in range(self.n_pairs):
            if self.pairwise_type[i] == 2:
                # Reduce excluded volume for really close LJ1210 contacts
                if self.pairwise_other_parameters[i][0] < 0.85:
                    if not (list(self.pairs[i]) in self.exclusions): 
                        self.exclusions.append(list(self.pairs[i]))
                        self.give_smaller_excluded_volume.append(list(self.pairs[i]))
            else: 
                # Always exclude pairs that aren't LJ1210.
                if list(self.pairs[i]) not in self.exclusions: 
                    self.exclusions.append(list(self.pairs[i]))

        if 4 in self.pairwise_type:
            self.contact_type = "Gaussian"
        else:
            self.contact_type = "LJ1210"

        # List of array indices indicating which interactions have the associated model parameter.
        self.n_model_param = len(self.model_param_values)
        self.model_param_interactions = [ (np.where(self.pairwise_param_assignment == p))[0] for p in range(self.n_model_param) ]

    def _determine_native_pairs(self):
        """ Determine which pairs are native pairs """

        # TODO(alex):
        # - Refactor/clean up this function.

        self.native_pairs_ndx = "[ native_contacts ]\n"
        self.native_pairs = []
        self.native_pairs_indices = []
        self.Qref = np.zeros((self.n_residues,self.n_residues))
        if self.n_native_pairs is None:
            if self.verbose:
                print "  Considering all (nonredundant) pairs to be native pairs"
            for i in range(self.n_pairs):
                if (list(self.pairs[i,:]) not in self.native_pairs) and (self.pairwise_type[i] not in SKIP_INTERACTIONS):
                    self.native_pairs_indices.append(i)
                    self.native_pairs.append(list(self.pairs[i,:]))
                    self.native_pairs_ndx += "%4d %4d\n" % (self.pairs[i,0],self.pairs[i,1])
                    if self.bead_repr == "CA":
                        self.Qref[self.pairs[i,0]-1,self.pairs[i,1]-1] = 1 
            self.n_native_pairs = len(self.native_pairs)
        else:
            if self.verbose:
                print "  Considering the first %d unique pairs to be native pairs" % self.n_native_pairs
            if self.n_native_pairs > self.n_pairs:
                raise IOError("%d n_native_pairs specified greater than %d n_pairs." % (self.n_native_pairs, self.n_pairs))

            # Takes first unique n_native_pairs as native.
            count = 0
            for i in range(self.n_pairs):
                if (list(self.pairs[i,:]) not in self.native_pairs) and (self.pairwise_type[i] not in SKIP_INTERACTIONS):
                    self.native_pairs_indices.append(i)
                    self.native_pairs.append(list(self.pairs[i,:]))
                    self.native_pairs_ndx += "%4d %4d\n" % (self.pairs[i,0],self.pairs[i,1])
                    count += 1
                    if self.bead_repr == "CA":
                        self.Qref[self.pairs[i,0]-1,self.pairs[i,1]-1] = 1 
                if count == self.n_native_pairs:
                    break
             
        if self.bead_repr == "CACB": 
            self.Qref = np.zeros((max(self.pairs.ravel())+1,max(self.pairs.ravel())+1))
            for i in range(self.n_pairs):
                self.Qref[self.pairs[i,0]-1,self.pairs[i,1]-1] = 1 
        
        self.native_pairs_indices = np.array(self.native_pairs_indices)

        # Estimate the stability of the native state
        self.native_stability = float(self.n_native_pairs)
        #self.native_stability = 0.
        #for i in range(self.n_native_pairs):
        #    idx = self.native_pairs_indices[i]
        #    self.native_stability += (self.pair_eps[idx]*self.pair_V[idx](np.array([self.pairwise_other_parameters[idx][0]])))[0]

    def _determine_tabled_interactions(self):
        """Determine which interactions need to use a tableb_##.xvg file"""
        int_type = self.pairwise_type
        if self.using_sbm_gmx:
            dont_table = [1,2,4,8,9,10]
        else:
            dont_table = [1,2]

        flag = np.zeros(len(int_type))
        for i in range(len(dont_table)):
            flag += (int_type == dont_table[i]).astype(int)

        self.tabled_interactions = np.zeros(self.n_pairs,float)
        for tbl_indx in (np.where(flag != 1))[0]:
            self.tabled_interactions[tbl_indx] = 1
        self.tabled_pairs = np.where(self.tabled_interactions == 1)[0]
        self.n_tables = len(self.tabled_pairs)

    def _set_nonbonded_interactions(self):
        """ Set all interaction functions """

        # File to save interaction parameters for pairwise potentials. 
        pair_eps = []
        self.pair_V = []
        self.pair_dV = []
        self.pairwise_param_file_string = "#   i   j   param int_type  other_params\n"
        for i in range(self.n_pairs):
            i_idx = self.pairs[i][0]
            j_idx = self.pairs[i][1]
            model_param = self.pairwise_param_assignment[i]
            int_type = self.pairwise_type[i]
            other_param_string = ""
            for p in range(len(self.pairwise_other_parameters[i])):
                other_param_string += " %10.5f " % self.pairwise_other_parameters[i][p] 
            self.pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (i_idx,j_idx,model_param,int_type,other_param_string)

            # Wrap the pairwise potentials so that only distance needs to be input.
            self.pair_V.append(pairwise.wrap_pairwise(pairwise.get_pair_potential(int_type),\
                                            *self.pairwise_other_parameters[i]))
            self.pair_dV.append(pairwise.wrap_pairwise(pairwise.get_pair_potential_deriv(int_type),\
                                            *self.pairwise_other_parameters[i]))

            # Assign pairwise interaction strength from model parameters
            pair_eps.append(self.model_param_values[self.pairwise_param_assignment[i]])

        self.pair_eps = np.array(pair_eps)

        # File to save model parameters
        self.model_param_file_string = "# model parameters\n"
        for i in range(self.n_model_param):
            self.model_param_file_string += "%10.5f\n" % self.model_param_values[i]

        if self.using_sbm_gmx:
            self._set_smog_pairs()
        else:
            self.smog_pair_indxs = []
        
        self._generate_interaction_tables()
        self._generate_topology()

    def _set_smog_pairs(self):
        """Create smog   """
        # Get unique pairs.
        smog_int_types = [4,8,9,10]
        self.smog_pairs = []
        self.smog_pair_indxs = []
        for i in range(self.n_pairs):
            p = self.pairs[i,:]
            if self.pairwise_type[i] in smog_int_types:
                if self.pairwise_type[i] == 8:
                    not_smog_pair = False
                    for j in range(i + 1,self.n_pairs):
                        p2 = self.pairs[j,:]
                        if np.all(p == p2) and (self.pairwise_type[j] == 5):
                            not_smog_pair = True
                            break
                    if not_smog_pair:
                        continue
                    else:
                        self.smog_pair_indxs.append(i)
                        if list(p) not in self.smog_pairs:
                            self.smog_pairs.append(list(p))
                else:
                    self.smog_pair_indxs.append(i)
                    if list(p) not in self.smog_pairs:
                        self.smog_pairs.append(list(p))

        self.smog_type = []
        self.smog_strings = []
        for i in range(len(self.smog_pairs)):
            pair = self.smog_pairs[i]
            temp = {}
            # Collect all interactions between given pair
            for j in range(self.n_pairs):
                if list(self.pairs[j]) == pair:
                    int_type = self.pairwise_type[j]  
                    if (int_type == 10) and (int_type in temp):
                        temp[-int_type] = j
                    else:
                        temp[int_type] = j

            # Determine their smog counterpart
            if 4 in temp:
                if 8 in temp: # Gaussian+LJ12 
                    ftype = 6
                    #r0, width = self.pairwise_other_parameters[temp[4]]
                    eps = self.pair_eps[temp[4]]
                    rNC, r0, width = self.pairwise_other_parameters[temp[8]]
                    smog_other = "%5d  %10.5e %10.5e %10.5e %10.5e" % (ftype,eps,r0,width,rNC**12)
                else: # Bare Gaussian
                    ftype = 5
                    r0, width = self.pairwise_other_parameters[temp[4]]
                    eps = self.pair_eps[temp[4]]
                    smog_other = "%5d  %10.e %10.5e %10.5e" % (ftype,eps,r0,width)
            else: # 2Gaussians+LJ12
                ftype = 7
                rNC, r0, width0, r1, width1 = self.pairwise_other_parameters[temp[9]]
                eps0 = self.pair_eps[temp[10]]
                #eps1 = self.pair_eps[temp[-10]]    # SBM only supports double well of equal depths.
                smog_other = "%5d  %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n" % (ftype,eps0,r0,width0,r1,width1,rNC**2)

            self.smog_type.append(ftype)
            self.smog_strings.append(smog_other)

    def _generate_interaction_tables(self):
        """ Generates tables of user-defined potentials 

        """
        # TODO(alex):
        #    - Sum all the tabled interactions between given pair? To reduce the # of table files.

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
            table[10:,1] = self.pair_V[pair_indx](r[10:]) 
            table[10:,2] = -1.*self.pair_dV[pair_indx](r[10:]) 
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
            pair = self.pairs[pair_indx]
            eps = self.pair_eps[pair_indx]
            tabled_string += "%6d %6d%2d%18d%18.9e\n" %  \
                      (pair[0],pair[1],9,i+1,eps)
        return tabled_string

    def _get_pairs_string(self):
        """ Get the [ pairs ] string """

        pairs_string = " [ pairs ]\n"
        pairs_string += " ; %5s  %5s %4s    %8s   %8s\n" % ("i","j","type","c10","c12")
        for i in range(self.n_pairs):
            if i not in self.smog_pair_indxs:
                res_a = self.pairs[i][0]
                res_b = self.pairs[i][1]
                r0 = self.pairwise_other_parameters[i][0]
                eps = self.pair_eps[i]
                if self.pairwise_type[i] == 1:     # LJ12
                    c12 = eps*(r0**12)
                    c10 = 0
                    pairs_string += "%6d %6d%2d%18.9e%18.9e\n" % (res_a,res_b,1,c10,c12)
                elif self.pairwise_type[i] == 2:   # LJ1210
                    c12 = eps*5.0*(r0**12)
                    c10 = eps*6.0*(r0**10)
                    pairs_string += "%6d %6d%2d%18.9e%18.9e\n" % (res_a,res_b,1,c10,c12)
                else:
                    pass

        # Give a smaller excluded volume for LJ1210 pairs that are closer.
        for i in range(len(self.give_smaller_excluded_volume)):
            res_a = self.give_smaller_excluded_volume[i][0]
            res_b = self.give_smaller_excluded_volume[i][1]
            c12 = (0.2**12)
            c10 = 0
            pairs_string += "%6d %6d%2d%18.9e%18.9e\n" % (res_a,res_b,1,c10,c12)
        
        if self.using_sbm_gmx:
            # Smog interactions
            pairs_string += " ; Smog guassian interactions\n"
            for i in range(len(self.smog_pairs)):
                res_a = self.smog_pairs[i][0]
                res_b = self.smog_pairs[i][1]
                pairs_string += "%5d%5d%s\n" % (res_a,res_b,self.smog_strings[i]) 

        return pairs_string

    def _get_pairs_string_smog(self):
        """ Get the [ pairs ] string for SBM Gromacs """

        # How to translate int_types to ftype?
        # Smog potential types:
        # Potential            ftype         int_type
        # Gaussian               5              4
        # Gaussian+LJ12          6              4+8
        # 2Gaussians+LJ12        7              9+10+10
        pairs_string = " [ pairs ]\n"
        pairs_string += " ; Using smog Guassian interactions\n"
        for i in range(len(self.smog_pairs)):
            res_a = self.smog_pairs[i][0]
            res_b = self.smog_pairs[i][1]
            pairs_string += "%5d%5d%s\n" % (res_a,res_b,self.smog_strings[i]) 

        return pairs_string

    def _get_exclusions_string(self):
        """ Get the [ exclusions ] string """
        if len(self.exclusions) > 0: 
            exclusions_string = " [ exclusions ]\n"
            exclusions_string += " ;  i    j \n"
            for i in range(len(self.exclusions)):
                res_a = self.exclusions[i][0]
                res_b = self.exclusions[i][1]
                exclusions_string += "%6d %6d\n" % (res_a,res_b)
        else:
            exclusions_string = ""

        return exclusions_string

    def _generate_topology(self):
        """ Return a structure-based topology file. Currently only for one molecule. """

        top_string =  " ; Structure-based  topology file for Gromacs:\n"
        top_string += " [ defaults ]\n"
        top_string += " ;nbfunc comb-rule gen-pairs\n"
        top_string += "      1           1 no\n\n"
        top_string += self.atomtypes_string
        top_string += " [ moleculetype ]\n"
        top_string += " ;name   nrexcl\n"
        top_string += " Macromolecule           3\n\n"

        top_string += self.atoms_string + "\n"
        top_string += self.bonds_string 
        top_string += self._get_tabled_string() + "\n"
        top_string += self.angles_string + "\n"
        top_string += self.dihedrals_string + "\n"
        top_string += self._get_pairs_string() + "\n"
        top_string += self._get_exclusions_string() + "\n"

        top_string += " [ system ]\n"
        top_string += " ; name\n"
        top_string += " Macromolecule\n\n"
        top_string += " [ molecules ]\n"
        top_string += " ; name molec \n"
        top_string += " Macromolecule 1\n\n"

        self.topology = top_string

    def save_simulation_files(self):
        """ Write all needed simulation files. """
        cwd = os.getcwd()
        relative_path = cwd.split("%s/" % self.path)[1]
        self.pairwise_params_file_location = "%s/pairwise_params" % relative_path
        self.model_params_file_location = "%s/model_params" % relative_path

        open("Native.pdb","w").write(self.cleanpdb)
        open("index.ndx","w").write(self.index_ndx)
        open("dihedrals.ndx","w").write(self.dihedrals_ndx)
        open("native_contacts.ndx","w").write(self.native_pairs_ndx)
        open("conf.gro","w").write(self.grofile)
        open("topol.top","w").write(self.topology)
        open("pairwise_params","w").write(self.pairwise_param_file_string)
        open("model_params","w").write(self.model_param_file_string)
        np.savetxt("Qref_cryst.dat",self.Qref,fmt="%1d",delimiter=" ")
        np.savetxt("pairs.dat",self.pairs,fmt="%4d",delimiter=" ")

        # Save needed table files
        np.savetxt("table.xvg",self.tablep,fmt="%16.15e",delimiter=" ")
        np.savetxt("tablep.xvg",self.tablep,fmt="%16.15e",delimiter=" ")
        for i in range(self.n_tables):
            np.savetxt(self.tablenames[i],self.tables[i],fmt="%16.15e",delimiter=" ")

    def update_model_param_values(self,new_model_param_values):
        """ If parameter changed sign, change the pairwise interaction type """
        # Switching between different interaction function types
        potential_type_switch = {2:3,3:2,4:5,5:4}
    
        # Loop over fitting_params only 
        for i in range(self.n_fitting_params):
            p_idx = self.fitting_params[i]
            p_pairs = self.model_param_interactions[p_idx]
            for n in range(len(p_pairs)):
                if new_model_param_values[i] < 0.:
                    # If parameter changes sign then the pairwise_type is flipped.
                    if self.pairwise_type[p_pairs[n]] == 1:
                        self.pairwise_type[p_pairs[n]] = 2
                    else:
                        self.pairwise_type[p_pairs[n]] = potential_type_switch[self.pairwise_type[p_pairs[n]]]
                else:
                    pass
                # Model parameters are always positive
                self.model_param_values[p_idx] = abs(new_model_param_values[i])   

        # Refresh everything that depends on model parameters
        self._determine_tabled_interactions()
        self._set_nonbonded_interactions()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build a SMOG model.')
    parser.add_argument('--name', type=str, required=True, help='pdb')
    parser.add_argument('--pairs', type=str, required=True, help='pairs')
    args = parser.parse_args() 

    name = args.name
    pdb = "%s.pdb" % name
    pairsfile = args.pairs

    if not os.path.exists(pairsfile):
        print "ERROR! file does not exists: ",pairsfile
        raise SystemExit
    else:
        pairs = np.loadtxt(pairsfile,dtype=int)

    model = CoarseGrainedModel(pdb=pdb,pairs=pairs,defaults=True,bead_repr="CACB")
    #model.save_simulation_files()
