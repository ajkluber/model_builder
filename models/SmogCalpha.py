''' SmogCalpha

Description:

    Generates topology and grofile needed to run Smog-style C-alpha Go-model. The
style of the input file formats emulates the SMOG server (see reference (1)).

Example Usage:
    See project_tools/examples

References:
'''

import numpy as np
import subprocess as sb
import os
import time

import bonded_potentials as bond
import pairwise_potentials as pairwise
import pdb_parser


class SmogCalpha(object):
    ''' This class creates a topology and grofile '''
    ## To Do:
    def __init__(self,**kwargs):

        self.path = os.getcwd()
        ## Set any keyword argument given as an attribute. Assumes it has what it needs.
        for key in kwargs.iterkeys():
            #print key.lower(),kwargs[key]
            setattr(self,key.lower(),kwargs[key])
        need_to_define = ["model_code","beadmodel","epsilon_bar",
                          "fitting_data","fitting_solver","fitting_allowswitch",
                          "disulfides","fitting_params"]
        for thing in need_to_define:
            if not hasattr(self,thing):
                setattr(self,thing,None)
        if not hasattr(self,"fitting_includes"):
            setattr(self,"fitting_includes",[self.pdb])
        if not hasattr(self,"defaults"):
            self.defaults = False

        if not os.path.exists(self.pdb):
            print "ERROR! The inputted pdb: %s does not exist" % pdb
            print " Exiting."
            raise SystemExit


        self.beadmodel = "CA"
        self.backbone_params = ["Kb","Ka","Kd"]
        self.backbone_param_vals = {"Kb":20000.,"Ka":400.,"Kd":1}

        self.error = 0
        self.initial_T_array = None
        self.name = self.pdb.split(".pdb")[0]
        self.subdir = self.name

        ## Set backbone parameters
        self._set_bonded_interactions()
        self.n_contacts = len(self.contacts)
        self.Qref = np.zeros((self.n_residues,self.n_residues))
        for pair in self.contacts:
            self.Qref[pair[0]-1,pair[1]-1] = 1 

        ## Check disulfide separation and remove from contacts list.
        self._check_disulfides() 

        ## Generate grofile, topology file, and necessary table files.
        self._check_contact_opts()
        self._get_interaction_tables()
        self._generate_index_ndx()
        self.generate_grofile()
        self.generate_topology()

    def get_model_info_string(self):
        ''' The string representation of all the model info.'''
        model_info_string = "[ Path ]\n"
        model_info_string += "%s\n" % self.path
        model_info_string += "[ PDB ]\n"
        model_info_string += "%s\n" % self.pdb
        model_info_string += "[ Subdir ]\n"
        model_info_string += "%s\n" % self.subdir
        model_info_string += "[ Iteration ]\n"
        model_info_string += "%s\n" % self.iteration
        model_info_string += "[ Model_Code ]\n" 
        model_info_string += "%s\n" % self.model_code
        model_info_string += "[ Bead_Model ]\n" 
        model_info_string += "%s\n" % self.beadmodel
        model_info_string += "[ Backbone_params ]\n" 
        model_info_string += "  %5s       %5s       %5s\n" % ("Kb","Ka","Kd")
        model_info_string += "[ Backbone_param_vals ]\n" 
        model_info_string += "%10.2f%10.2f%10.2f\n" % \
            (self.backbone_param_vals["Kb"],self.backbone_param_vals["Ka"],self.backbone_param_vals["Kd"])
        model_info_string += "[ Disulfides ]\n"
        if self.disulfides == None:
            model_info_string += "%s\n" % None
        else:
            temp = ''
            for x in self.disulfides:
                temp += " %d " % x 
            model_info_string += "%s\n" % temp
        model_info_string += "[ Epsilon_Bar ]\n"
        model_info_string += "%s\n" % str(self.epsilon_bar)
        model_info_string += "[ Nonnative ]\n"
        model_info_string += "%s\n" % str(None)
        model_info_string += "[ Pairwise_Params_File ]\n"
        model_info_string += "%s\n" % "None"
        model_info_string += "[ Model_Params_File ]\n"
        model_info_string += "%s\n" % "None"
        model_info_string += "[ Fitting_Data ]\n"
        model_info_string += "%s\n" % self.fitting_data
        model_info_string += "[ Fitting_Includes ]\n"
        for dir in self.fitting_includes:
            model_info_string += "%s" % str(dir)
        model_info_string += "\n"
        model_info_string += "[ Fitting_Solver ]\n"
        model_info_string += "%s\n" % self.fitting_solver
        model_info_string += "[ Fitting_AllowSwitch ]\n"
        model_info_string += "%s\n" % self.fitting_allowswitch
        model_info_string += "[ Fitting_Params ]\n"
        model_info_string += "%s\n" % self.fitting_params
        model_info_string += "[ Reference ]\n" 
        model_info_string += "None\n" 
        return model_info_string

    def append_log(self,string):
        now = time.localtime()
        now_string = "%s:%s:%s:%s:%s" % (now.tm_year,now.tm_mon,now.tm_mday,now.tm_hour,now.tm_min)
        logfile = open('%s/%s/%s.log' % (path,self.subdir,self.subdir),'a').write("%s %s\n" % (now_string,string))

    def _check_contact_opts(self):
        ''' Set default pairwise interaction terms '''
        ## All non-native pairs
        self.nonnative_pairs = []
        for i in range(self.n_residues):
            for j in range(i+4,self.n_residues):
                self.nonnative_pairs.append([i+1,j+1])

        self.nearnative_pairs = list(self.nonnative_pairs)  ## To Do: List of near-native contacts
        ## Remove native pairs
        for n in range(self.n_contacts):
            self.nonnative_pairs.pop(self.nonnative_pairs.index(list(self.contacts[n,:])))

        ## Grab structural distances.
        self.pairwise_distances = np.zeros(len(self.contacts),float)
        for i in range(len(self.contacts)):
            i_idx = self.contacts[i][0]
            j_idx = self.contacts[i][1]
            self.pairwise_distances[i] = bond.distance(self.atom_coords,i_idx-1,j_idx-1)

        ## Set some defaults for pairwise interaction potential.
        if not hasattr(self,"epsilon_bar"):
            self.epsilon_bar = None

        ## If defaults is False then all these need to be defined.
        if self.defaults:
            print "  Using defaults: LJ1210 contacts, homogeneous contacts 1, each contact free param."
            self.model_param_values = np.ones(self.n_contacts,float)
            self.pairwise_type = 2*np.ones(self.n_contacts,float)
            self.pairwise_param_assignment = np.arange(self.n_contacts)     
            self.pairwise_other_parameters = [ [self.pairwise_distances[x]] for x in range(self.n_contacts) ]
        else:
            needed = [hasattr(self,"model_param_values"),hasattr(self,"pairwise_type"), \
                    hasattr(self,"pairwise_param_assignment"),hasattr(self,"pairwise_other_parameters")]
            if not all(needed):
                print "ERROR! If not using  "

        self.n_model_param = len(self.model_param_values)
            
        self.tabled_interactions = np.zeros(self.n_contacts,float)
        for rep_indx in (np.where(self.pairwise_type != 2))[0]:
            self.tabled_interactions[rep_indx] = 1
        self.tabled_pairs = np.where(self.tabled_interactions == 1)[0]
        self.n_tables = len(self.tabled_pairs)

        ## List of array indices indicating which interactions have the associated model parameter.
        self.model_param_interactions = [ np.where(self.pairwise_param_assignment == p) for p in range(self.n_model_param) ]
        self.pairwise_strengths = np.array([  self.model_param_values[x] for x in self.pairwise_param_assignment ])

        ## Wrap the pairwise contact potentials so that only distance needs to be input.
        self.pairwise_potentials = [ pairwise.wrap_pairwise(pairwise.get_pair_potential(self.pairwise_type[x]),\
                                                *self.pairwise_other_parameters[x]) for x in range(self.n_contacts) ]

        self.pairwise_potentials_deriv = [ pairwise.wrap_pairwise(pairwise.get_pair_potential_deriv(self.pairwise_type[x]),\
                                                *self.pairwise_other_parameters[x]) for x in range(self.n_contacts) ]
        ## File to save model parameters
        self.model_param_file = "# model parameters\n"
        for i in range(self.n_model_param):
            self.model_param_file += "%10.5f\n" % self.model_param_values[i]

        ## File to save interaction parameters for pairwise potentials. 
        self.pairwise_param_file = "#   i   j   param int_type  other_params\n"
        for i in range(self.n_contacts):
            i_idx = self.contacts[i][0]
            j_idx = self.contacts[i][1]
            model_param = self.pairwise_param_assignment[i]
            int_type = self.pairwise_type[i]
            other_param_string = ""
            for p in range(len(self.pairwise_other_parameters[i])):
                other_param_string += " %10.5f " % self.pairwise_other_parameters[i][p] 
            
            self.pairwise_param_file += "%5d%5d%5d%5d%s\n" % (i_idx,j_idx,model_param,int_type,other_param_string)

    def update_pairwise_interactions(self):
        ''' Update pairwise interactions 

            Update pairwise interactions in case the interaction type or
        interaction strengths has changed.

        NOT DONE

        '''

        ## List of array indices indicating which interactions have the associated model parameter.
        self.model_param_interactions = [ np.where(self.pairwise_param_assignment == p) for p in range(self.n_model_param) ]
        self.pairwise_strengths = np.array([  self.model_param_values[x] for x in self.pairwise_param_assignment ])

        ## Wrap the pairwise contact potentials so that only distance needs to be input.
        self.pairwise_potentials = [ pairwise.wrap_pairwise(pairwise.get_pair_potential(self.pairwise_type[x]),\
                                                *self.pairwise_other_parameters[x]) for x in range(self.n_contacts) ]
        self.pairwise_potentials_deriv = [ pairwise.wrap_pairwise(pairwise.get_pair_potential_deriv(self.pairwise_type[x]),\
                                                *self.pairwise_other_parameters[x]) for x in range(self.n_contacts) ]
        ## File to save model parameters
        self.model_param_file = "# model parameters\n"
        for i in range(self.n_model_param):
            self.model_param_file += "%10.5f\n" % self.model_param_values[i]

        ## File to save interaction parameters for pairwise potentials. 
        self.pairwise_param_file = "#   i   j   param int_type  other_params\n"
        for i in range(self.n_contacts):
            i_idx = self.contacts[i][0]
            j_idx = self.contacts[i][1]
            model_param = self.pairwise_param_assignment[i]
            int_type = self.pairwise_type[i]
            other_param_string = ""
            for p in range(len(self.pairwise_other_parameters[i])):
                other_param_string += " %10.5f " % self.pairwise_other_parameters[i][p]

            self.pairwise_param_file += "%5d%5d%5d%5d%s\n" % (i_idx,j_idx,model_param,int_type,other_param_string)


    def _check_disulfides(self):
        ''' Check that specified disulfides are between cysteine and that 
            the corresonding pairs are within 0.8 nm.

        To Do: 
            - Add disulfides to the list of atoms to be bonded.
        '''
        coords = self.atom_coords
        residues = self.atom_residues
        if self.disulfides != None:
            print "  Checking if disulfides are reasonable."
            for i in range(len(self.disulfides[::2])):
                i_idx = self.disulfides[2*i]
                j_idx = self.disulfides[2*i + 1]
                dist = bond.distance(coords,i_idx-1,j_idx-1)
                theta1 = bond.angle(coords,i_idx-2,i_idx-1,j_idx-1)
                theta2 = bond.angle(coords,i_idx-1,j_idx-1,j_idx-2)
                phi = bond.dihedral(coords,i_idx-2,i_idx-1,j_idx-1,j_idx-2)
                if (residues[i_idx-1] != "CYS") or (residues[j_idx-1] != "CYS"):
                    print "WARNING! Specifying disulfide without cysteines: %s  %s" % \
                            (residues[i_idx-1]+str(i_idx), residues[j_idx-1]+str(j_idx))
                if dist > 0.8:
                    print "WARNING! Specifying disulfide with separation greater than 0.8 nm."
                else:
                    print "   %s %s separated by %.4f nm, Good." % \
                    (residues[i_idx-1]+str(i_idx), residues[j_idx-1]+str(j_idx),dist)
                    ## Remove disulfide pair from self.contacts if it is there.
                    new_conts = []
                    for pair in self.contacts:
                        if (pair[0] == i_idx) and (pair[1] == j_idx):
                            continue
                        else:
                            new_conts.append(pair)
                    self.contacts = np.array(new_conts)

                    ## Set cysteine bond distance, angles, and dihedral.
                    self.bond_indices.append([i_idx,j_idx])
                    self.bond_min.append(dist)
                    self.bond_strengths.append(self.backbone_param_vals["Kb"])

                    self.angle_indices.append([i_idx-1,i_idx,j_idx])
                    self.angle_indices.append([i_idx,j_idx,j_idx-1])
                    self.angle_min.append(theta1)
                    self.angle_min.append(theta2)
                    self.angle_strengths.append(self.backbone_param_vals["Ka"])
                    self.angle_strengths.append(self.backbone_param_vals["Ka"])

                    self.dihedral_indices.append([i_idx-1,i_idx,j_idx,j_idx-1])
                    self.dihedral_min.append(phi)
                    self.dihedral_strengths.append(self.backbone_param_vals["Kd"])

                if self.Qref[i_idx-1][j_idx-1] == 1:
                    #print "    Subtracting 1 from n_contacts for ", i_idx,j_idx, " disulfide"
                    self.n_contacts -= 1
        else:
            print "  No disulfides to check."

        self.contacts_ndx = "[ contacts ]\n"
        for i in range(self.n_contacts):
            self.contacts_ndx += "%4d %4d\n" % (self.contacts[i][0],self.contacts[i][1])

    def _set_bonded_interactions(self):
        ''' Extract info from the Native.pdb for making index and top file '''
        ## Grab coordinates from the pdb file.
        self.cleanpdb_full, self.cleanpdb_full_noH, self.cleanpdb = pdb_parser.clean(self.pdb)
        coords, indices, atoms, residues = pdb_parser.get_coords_atoms_residues(self.cleanpdb)

        self.atom_indices = indices
        self.atom_types = atoms
        self.atom_residues = residues
        self.atom_coords = coords

        self.n_residues = len(residues)
        self.n_atoms = len(atoms)

        ## Set bonded force field terms quantities.
        self.bond_indices = [[indices[i],indices[i+1]] for i in range(self.n_atoms-1)]
        self.angle_indices = [[indices[i],indices[i+1],indices[i+2]] for i in range(self.n_atoms-2)]
        self.dihedral_indices = [[indices[i],indices[i+1],indices[i+2],indices[i+3]] for i in range(self.n_atoms-3)]

        self.bond_min = [ bond.distance(coords,i_idx-1,j_idx-1) for i_idx,j_idx in self.bond_indices ]
        self.angle_min = [ bond.angle(coords,i_idx-1,j_idx-1,k_idx-1) for i_idx,j_idx,k_idx in self.angle_indices ]
        self.dihedral_min = [ bond.dihedral(coords,i_idx-1,j_idx-1,k_idx-1,l_idx-1) for i_idx,j_idx,k_idx,l_idx in self.dihedral_indices ]

        if not hasattr(self,"bond_strengths"):
            self.bond_strengths = [ self.backbone_param_vals["Kb"] for i in range(len(self.bond_min)) ]
        if not hasattr(self,"angle_strengths"):
            self.angle_strengths = [ self.backbone_param_vals["Ka"] for i in range(len(self.angle_min)) ]
        if not hasattr(self,"dihedral_strengths"):
            self.dihedral_strengths = [ self.backbone_param_vals["Kd"] for i in range(len(self.dihedral_min)) ]

    def _get_interaction_tables(self):
        ''' Generates tables of user-defined potentials '''

        self.tablep = self._get_LJ1210_table()
        self.LJtable = self.tablep
        r = np.arange(0.002,100.0,0.002)
        self.tables = []
        self.tablenames = []
        for i in range(self.n_tables):
            pair_indx = self.tabled_pairs[i]
            table_name = "table_b%d.xvg" % (i+1)
            self.tablenames.append(table_name)

            table = np.zeros((len(r)+1,3),float)
            table[1:,0] = r 
            table[1:,1] = self.pairwise_potentials[pair_indx](r) 
            table[1:,2] = -1.*self.pairwise_potentials_deriv[pair_indx](r) 
            self.tables.append(table)

    def _get_LJ1210_table(self):
        ''' LJ12-10 interaction potential ''' 
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

    def _generate_index_ndx(self):
        ''' Generates index file for gromacs analysis utilities. '''
        ca_string = ''
        i = 1
        
        for indx in self.atom_indices: 
            if (i % 15) == 0:
                ca_string += '%4d \n' % indx
            else:
                ca_string += '%4d ' % indx
            i += 1
        ca_string += '\n'
        headings = ["system","Protein","Protein-H","C-alpha",\
                    "Backbone","MainChain","MainChain+Cb","MainChain+H"]
        indexstring = ""
        for heading in headings:
            indexstring += "[ "+heading+" ]\n"
            indexstring += ca_string
        indexstring += '[ SideChain ]\n\n'
        indexstring += '[ SideChain-H ]\n\n'
        self.index_ndx = indexstring

    def _get_atoms_string(self):
        ''' Generate the [ atoms ] string.'''
        atoms_string = " [ atoms ]\n"
        atoms_string += " ;nr  type  resnr residue atom  cgnr charge  mass\n"
        for j in range(len(self.atom_indices)):
            atomnum = self.atom_indices[j]
            resnum = self.atom_indices[j]
            resid = self.atom_residues[j]
            atoms_string += " %5d%4s%8d%5s%4s%8d%8.3f%8.3f\n" % \
                        (atomnum,"CA",resnum,resid,"CA",atomnum,0.0,1.0)
        return atoms_string

    def _get_bonds_string(self):
        ''' Generate the [ bonds ] string.'''
        bonds_string = " [ bonds ]\n"
        bonds_string += " ; ai aj func r0(nm) Kb\n"
        for j in range(len(self.bond_min)):
            i_idx = self.bond_indices[j][0]
            j_idx = self.bond_indices[j][1]
            dist = self.bond_min[j]
            kb = self.bond_strengths[j]
            bonds_string += "%6d %6d%2d%18.9e%18.9e\n" %  \
                          (i_idx,j_idx,1,dist,kb)
        return bonds_string

    def _get_tabled_string(self):
        ''' Generate the topology files to specify table interactions. '''
        ## Add special nonbonded table interactions. 
        tabled_string = "; tabled interactions contacts below\n"
        for i in range(self.n_tables):
            pair = self.contacts[self.tabled_pairs[i]]
            tabled_string += "%6d %6d%2d%18d%18.9e\n" %  \
                      (pair[0],pair[1],9,i+1,1.0)
        return tabled_string

    def _get_angles_string(self):
        ''' Generate the [ angles ] string.'''
        angles_string = " [ angles ]\n"
        angles_string += " ; ai  aj  ak  func  th0(deg)   Ka\n"
        for j in range(len(self.angle_min)):
            i_idx = self.angle_indices[j][0]
            j_idx = self.angle_indices[j][1]
            k_idx = self.angle_indices[j][2]
            theta = self.angle_min[j]
            ka = self.angle_strengths[j]
            angles_string += "%6d %6d %6d%2d%18.9e%18.9e\n" %  \
                          (i_idx,j_idx,k_idx,1,theta,ka)
        return angles_string

    def _get_dihedrals_string(self):
        ''' Generate the [ dihedrals ] string.'''
        dihedrals_string = " [ dihedrals ]\n"
        dihedrals_string += " ; ai  aj  ak al  func  phi0(deg)   Kd mult\n"
        dihedrals_ndx = '[ dihedrals ]\n'
        for j in range(len(self.dihedral_min)):
            i_idx = self.dihedral_indices[j][0]
            j_idx = self.dihedral_indices[j][1]
            k_idx = self.dihedral_indices[j][2]
            l_idx = self.dihedral_indices[j][3]
            phi = self.dihedral_min[j]
            kd = self.dihedral_strengths[j]
            dihedrals_string += "%6d %6d %6d %6d%2d%18.9e%18.9e%2d\n" %  \
                          (i_idx,j_idx,k_idx,l_idx,1,phi,kd,1)
            dihedrals_string += "%6d %6d %6d %6d%2d%18.9e%18.9e%2d\n" %  \
                          (i_idx,j_idx,k_idx,l_idx,1,3.*phi,kd/2.,3)
            dihedrals_ndx += '%4d %4d %4d %4d\n' % \
                                (i_idx,j_idx,k_idx,l_idx)
        self.dihedrals_ndx = dihedrals_ndx
        return dihedrals_string

    def _get_pairs_string(self):
        ''' Get the [ pairs ] string '''

        if self.epsilon_bar != None:
            print "  Scaling attractive contacts such that epsilon_bar =", self.epsilon_bar
            ## TO DO:
            ## - epsilon bar for attractive contacts
            #stability = sum(self.contact_epsilons[self.LJtype == 1]) - sum(
            #avg_eps = np.mean(self.contact_epsilons[attractive])
            #self.contact_epsilons[attractive] = /avg_eps

        pairs_string = " [ pairs ]\n"
        #pairs_string += " ; i    j      type  c10  \n"
        pairs_string += " ; %5s  %5s %4s    %8s   %8s\n" % ("i","j","type","c10","c12")
        beadbead_string = ""
        for i in range(self.n_contacts):
            res_a = self.contacts[i][0]
            res_b = self.contacts[i][1]
            resid_a = self.atom_residues[res_a-1]
            resid_b = self.atom_residues[res_b-1]
            sig_ab = self.pairwise_distances[i]
            eps_ab = self.pairwise_strengths[i]

            if self.pairwise_type[i] == 2:   ## LJ1210
                c12 = eps_ab*5.0*(sig_ab**12)
                c10 = eps_ab*6.0*(sig_ab**10)

                pairs_string += "%6d %6d%2d%18.9e%18.9e\n" % (res_a,res_b,1,c10,c12)
                beadbead_string += "%5d%5d%8s%8s%5d%18.9e%18.9e%18d\n" % \
                                (res_a,res_b,resid_a,resid_b,i+1,sig_ab,eps_ab,self.pairwise_type[i])

            #elif self.pairwise_type[i] == 4: ## Gaussian
            #    width_ab = self.contact_widths[i] 
            #    ncwall_ab = self.noncontact_wall[i]**12
            #    pairs_string += "%6d %6d%2d%18.9e%18.9e%18.9e%18.9e\n" % \
            #                    (res_a,res_b,6,eps_ab,sig_ab,width_ab,ncwall_ab)
            #    beadbead_string += "%5d%5d%8s%8s%5d%18.9e%18.9e%18d%18.9e%18.9e\n" % \
            #        (res_a,res_b,resid_a,resid_b,i+1,sig_ab,eps_ab,1,width_ab,ncwall_ab)
            #elif self.pairwise_type[i] == 4: ## Gaussian
            #    width_ab = self.contact_widths[i] 
            #    ncwall_ab = self.noncontact_wall[i]**12
            #    pairs_string += "%6d %6d%2d%18.9e%18.9e%18.9e%18.9e\n" % \
            #                    (res_a,res_b,6,eps_ab,sig_ab,width_ab,ncwall_ab)
            #    beadbead_string += "%5d%5d%8s%8s%5d%18.9e%18.9e%18d%18.9e%18.9e\n" % \
            #        (res_a,res_b,resid_a,resid_b,i+1,sig_ab,eps_ab,1,width_ab,ncwall_ab)
            else:
                continue
                #print "ERROR! unrecognized contact option: ", self.pairwise_type[i]
                #print "Exiting"
                #raise SystemExit

        self.beadbead = beadbead_string 
        return pairs_string

    def _get_exclusions_string(self):
        ''' Get the [ exclusions ] string '''
        exclusions_string = ""
        #exclusions_string = " [ exclusions ]\n"
        #exclusions_string += " ; ai aj \n"
        ## Removing exlusions for pairwise interactions so that hard-sphere 
        ## radius of beads is always respected.
        #for i in range(len(self.contacts)):
        #    res_a = self.contacts[i][0]
        #    res_b = self.contacts[i][1]
        #    exclusions_string += "%6d %6d\n" % (res_a,res_b)
        if self.disulfides != None:
            exclusions_string = " [ exclusions ]\n"
            exclusions_string += " ; ai aj \n"
            for i in range(len(self.disulfides)/2):
                cys_a = self.disulfides[2*i]
                cys_b = self.disulfides[2*i + 1]
                exclusions_string += "%6d %6d\n" % (cys_a,cys_b)

        return exclusions_string

    def generate_grofile(self):
        ''' Get the .gro string '''
        gro_string = " Structure-based Gro file\n"
        gro_string += "%12d\n" % len(self.atom_types)
        for i in range(len(self.atom_types)):
            gro_string += "%5d%5s%4s%6d%8.3f%8.3f%8.3f\n" % \
                (self.atom_indices[i],self.atom_residues[i],self.atom_types[i],self.atom_indices[i],
                self.atom_coords[i][0],self.atom_coords[i][1],self.atom_coords[i][2])
        gro_string += "   %-25.16f%-25.16f%-25.16f" % (100.0,100.0,100.0)

        self.grofile = gro_string

    def generate_topology(self):
        ''' Return a structure-based topology file. Currently only for one molecule. '''

        top_string =  " ; Structure-based  topology file for Gromacs:\n"
        top_string += " [ defaults ]\n"
        top_string += " ;nbfunc comb-rule gen-pairs\n"
        top_string += "      1           1 no\n\n"
        top_string += " [ atomtypes ]\n"
        top_string += " ;name  mass     charge   ptype c10       c12\n"
        #top_string += " CA     1.000    0.000 A    0.000   0.167772160E-04\n\n"
        top_string += " CA     1.000    0.000 A    0.000   %10.9e\n\n" % (0.2**12) 
        top_string += " [ moleculetype ]\n"
        top_string += " ;name   nrexcl\n"
        top_string += " Macromolecule           3\n\n"

        top_string += self._get_atoms_string() + "\n"
        top_string += self._get_bonds_string() 
        top_string += self._get_tabled_string() + "\n"
        top_string += self._get_angles_string() + "\n"
        top_string += self._get_dihedrals_string() + "\n"
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
        ''' Write all needed simulation files. '''
        open("Native.pdb","w").write(self.cleanpdb)
        open("index.ndx","w").write(self.index_ndx)
        open("dihedrals.ndx","w").write(self.dihedrals_ndx)
        open("contacts.ndx","w").write(self.contacts_ndx)
        open("conf.gro","w").write(self.grofile)
        open("topol.top","w").write(self.topology)
        open("pairwise_params","w").write(self.pairwise_param_file)
        open("model_params","w").write(self.model_param_file)
        np.savetxt("Qref_cryst.dat",self.Qref,fmt="%1d",delimiter=" ")
        np.savetxt("contacts.dat",self.contacts,fmt="%4d",delimiter=" ")

        ## Save needed table files
        np.savetxt("table.xvg",self.tablep,fmt="%16.15e",delimiter=" ")
        for i in range(self.n_tables):
            np.savetxt(self.tablenames[i],self.tables[i],fmt="%16.15e",delimiter=" ")

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='Build a SMOG model.')
    parser.add_argument('--name', type=str, required=True, help='pdb')
    parser.add_argument('--contacts', type=str, required=True, help='contacts')
    args = parser.parse_args() 

    name = args.name
    pdb = "%s.pdb" % name
    contactsfile = args.contacts

    if not os.path.exists(contactsfile):
        print "ERROR! file does not exists: ",contactsfile
        raise SystemExit
    else:
        contacts = np.loadtxt(contactsfile,dtype=int)

    model = SmogCalpha(pdb=pdb,contacts=contacts)
    model.save_simulation_files()
