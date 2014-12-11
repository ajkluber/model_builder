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
    def __init__(self,**kwargs):
        self.path = os.getcwd()
        ## Set any keyword argument given as an attribute. Assumes it has what it needs.
        for key in kwargs.iterkeys():
            #print key.lower(),kwargs[key]
            setattr(self,key.lower(),kwargs[key])
        need_to_define = ["model_code","beadmodel","epsilon_bar",
                          "fitting_data","fitting_solver","fitting_allowswitch",
                          "disulfides","fitting_params","nonnative",
                          "pairwise_params_file_location","model_params_file_location"]
        for thing in need_to_define:
            if not hasattr(self,thing):
                setattr(self,thing,None)
        if not hasattr(self,"fitting_includes"):
            self.fitting_includes = [self.pdb]
        if not hasattr(self,"defaults"):
            self.defaults = False

        if not os.path.exists(self.pdb):
            print "ERROR! The inputted pdb: %s does not exist" % pdb
            print " Exiting."
            raise SystemExit

        self.beadmodel = "CA"
        self.backbone_params = ["Kb","Ka","Kd"]
        self.backbone_param_vals = {"Kb":20000.,"Ka":400.,"Kd":1}

        self.initial_T_array = None
        self.name = self.pdb.split(".pdb")[0]
        self.subdir = self.name
        self.exclusions = []

        ## Set backbone parameters
        self._set_bonded_interactions()
        self._generate_index_ndx()
        self._generate_grofile()
        self.n_contacts = len(self.contacts)
        self.Qref = np.zeros((self.n_residues,self.n_residues))
        for pair in self.contacts:
            self.Qref[pair[0]-1,pair[1]-1] = 1 

        ## Check disulfide separation and remove from contacts list.
        self._check_disulfides() 

        ## Generate topology file and table files.
        self._check_contact_opts()
        self._set_nonbonded_interactions()

    def get_model_info_string(self):
        ''' The string representation of all the model info.'''
        model_info_string = "[ Path ]\n"
        model_info_string += "%s\n" % self.path
        model_info_string += "[ PDB ]\n"
        model_info_string += "%s\n" % self.pdb
        model_info_string += "[ Subdir ]\n"
        model_info_string += "%s\n" % self.subdir
        model_info_string += "[ Iteration ]\n"
        model_info_string += "%d\n" % self.iteration
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
        model_info_string += "%s\n" % str(self.pairwise_params_file_location)
        model_info_string += "[ Model_Params_File ]\n"
        model_info_string += "%s\n" % str(self.model_params_file_location)
        model_info_string += "[ Contact_Type ]\n"
        model_info_string += "%s\n" % self.contact_type
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
        self.n_nonnative_contacts = len(self.nonnative_pairs)

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
                print "ERROR! If not using defualts then the following need to be inputted:"
                print "     model_param_values"
                print "     pairwise_type"
                print "     pairwise_param_assignment"
                print "     pairwise_other_parameters"
        

        ## Any contact that would be disrupted by an 
        ## excluded volume radius of 0.4nm
        for i in range(self.n_contacts):
            if self.pairwise_other_parameters[i][0] < 0.85:
                self.exclusions.append(list(self.contacts[i]))

        self.contact_type = "none"
        for i in self.pairwise_type:
            if i == 4:
                self.contact_type="Gaussian"
        
        if self.contact_type == "none":
            self.contact_type = "LJ1210"

        ## TO DO: - How to scale the model parameters to get constant stability?
        if self.epsilon_bar != None:
            pass

        ## List of array indices indicating which interactions have the associated model parameter.
        self.n_model_param = len(self.model_param_values)
        self.model_param_interactions = [ np.where(self.pairwise_param_assignment == p) for p in range(self.n_model_param) ]

    def _check_disulfides(self):
        ''' Check that specified disulfides are between cysteine and that 
            the corresonding pairs are within 0.8 nm. '''
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

                    self.exclusions.append([i_idx,j_idx])
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
            
    def _set_nonbonded_interactions(self):
        ''' Set all interaction functions '''

        ## Determine the number of tabled interactions. Need to table if not LJ1210 or LJ12
        flag = ((self.pairwise_type == 2).astype(int) + (self.pairwise_type == 1).astype(int))
        self.tabled_interactions = np.zeros(self.n_contacts,float)
        for rep_indx in (np.where(flag != 1))[0]:
            self.tabled_interactions[rep_indx] = 1
        self.tabled_pairs = np.where(self.tabled_interactions == 1)[0]
        self.n_tables = len(self.tabled_pairs)

        ## Assign pairwise interaction strength from model parameters
        self.pairwise_strengths = np.array([  self.model_param_values[x] for x in self.pairwise_param_assignment ])

        ## Wrap the pairwise contact potentials so that only distance needs to be input.
        self.pairwise_potentials = [ pairwise.wrap_pairwise(pairwise.get_pair_potential(self.pairwise_type[x]),\
                                                *self.pairwise_other_parameters[x]) for x in range(self.n_contacts) ]

        self.pairwise_potentials_deriv = [ pairwise.wrap_pairwise(pairwise.get_pair_potential_deriv(self.pairwise_type[x]),\
                                                *self.pairwise_other_parameters[x]) for x in range(self.n_contacts) ]

        ## File to save model parameters
        self.model_param_file_string = "# model parameters\n"
        for i in range(self.n_model_param):
            self.model_param_file_string += "%10.5f\n" % self.model_param_values[i]

        ## File to save interaction parameters for pairwise potentials. 
        self.pairwise_param_file_string = "#   i   j   param int_type  other_params\n"
        for i in range(self.n_contacts):
            i_idx = self.contacts[i][0]
            j_idx = self.contacts[i][1]
            model_param = self.pairwise_param_assignment[i]
            int_type = self.pairwise_type[i]
            other_param_string = ""
            for p in range(len(self.pairwise_other_parameters[i])):
                other_param_string += " %10.5f " % self.pairwise_other_parameters[i][p] 
            self.pairwise_param_file_string += "%5d%5d%5d%5d%s\n" % (i_idx,j_idx,model_param,int_type,other_param_string)

        self._generate_interaction_tables()
        self._generate_topology()

    def _generate_interaction_tables(self):
        ''' Generates tables of user-defined potentials '''

        self.tablep = self._get_LJ1210_table()
        self.LJtable = self.tablep
        r = np.arange(0,100.0,0.002)
        self.tables = []
        self.tablenames = []
        for i in range(self.n_tables):
            pair_indx = self.tabled_pairs[i]
            table_name = "table_b%d.xvg" % (i+1)
            self.tablenames.append(table_name)

            table = np.zeros((len(r),3),float)
            table[1:,0] = r[1:]
            table[10:,1] = self.pairwise_potentials[pair_indx](r[10:]) 
            table[10:,2] = -1.*self.pairwise_potentials_deriv[pair_indx](r[10:]) 
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

        pairs_string = " [ pairs ]\n"
        pairs_string += " ; %5s  %5s %4s    %8s   %8s\n" % ("i","j","type","c10","c12")
        for i in range(self.n_contacts):
            res_a = self.contacts[i][0]
            res_b = self.contacts[i][1]
            r0 = self.pairwise_distances[i]
            eps = self.pairwise_strengths[i]
            if self.pairwise_type[i] == 1:     ## LJ12
                c12 = eps*(r0**12)
                c10 = 0
                pairs_string += "%6d %6d%2d%18.9e%18.9e\n" % (res_a,res_b,1,c10,c12)
            elif self.pairwise_type[i] == 2:   ## LJ1210
                c12 = eps*5.0*(r0**12)
                c10 = eps*6.0*(r0**10)
                pairs_string += "%6d %6d%2d%18.9e%18.9e\n" % (res_a,res_b,1,c10,c12)
            else:
                pass

        ## Give a smaller hard-wall exclusion for contacts that are closer.
        for i in range(len(self.exclusions)):
            res_a = self.exclusions[i][0]
            res_b = self.exclusions[i][1]
            c12 = (0.2**12)
            c10 = 0
            pairs_string += "%6d %6d%2d%18.9e%18.9e\n" % (res_a,res_b,1,c10,c12)

        return pairs_string

    def _get_exclusions_string(self):
        ''' Get the [ exclusions ] string '''
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

    def _generate_grofile(self):
        ''' Get the .gro string '''
        gro_string = " Structure-based Gro file\n"
        gro_string += "%12d\n" % len(self.atom_types)
        for i in range(len(self.atom_types)):
            gro_string += "%5d%5s%4s%6d%8.3f%8.3f%8.3f\n" % \
                (self.atom_indices[i],self.atom_residues[i],self.atom_types[i],self.atom_indices[i],
                self.atom_coords[i][0],self.atom_coords[i][1],self.atom_coords[i][2])
        gro_string += "   %-25.16f%-25.16f%-25.16f" % (100.0,100.0,100.0)

        self.grofile = gro_string

    def _generate_topology(self):
        ''' Return a structure-based topology file. Currently only for one molecule. '''

        top_string =  " ; Structure-based  topology file for Gromacs:\n"
        top_string += " [ defaults ]\n"
        top_string += " ;nbfunc comb-rule gen-pairs\n"
        top_string += "      1           1 no\n\n"
        top_string += " [ atomtypes ]\n"
        top_string += " ;name  mass     charge   ptype c10       c12\n"
        top_string += " CA     1.000    0.000 A    0.000   %10.9e\n\n" % (0.4**12) ## 0.4nm excluded vol radius
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
        cwd = os.getcwd()
        self.pairwise_params_file_location = "%s/pairwise_params" % cwd
        self.model_params_file_location = "%s/model_params" % cwd

        open("Native.pdb","w").write(self.cleanpdb)
        open("index.ndx","w").write(self.index_ndx)
        open("dihedrals.ndx","w").write(self.dihedrals_ndx)
        open("contacts.ndx","w").write(self.contacts_ndx)
        open("conf.gro","w").write(self.grofile)
        open("topol.top","w").write(self.topology)
        open("pairwise_params","w").write(self.pairwise_param_file_string)
        open("model_params","w").write(self.model_param_file_string)
        np.savetxt("Qref_cryst.dat",self.Qref,fmt="%1d",delimiter=" ")
        np.savetxt("contacts.dat",self.contacts,fmt="%4d",delimiter=" ")

        ## Save needed table files
        np.savetxt("tablep.xvg",self.tablep,fmt="%16.15e",delimiter=" ")
        for i in range(self.n_tables):
            np.savetxt(self.tablenames[i],self.tables[i],fmt="%16.15e",delimiter=" ")

    def update_model_param_values(self,new_model_param_values):
        ''' If parameter changed sign, change the pairwise interaction type '''

        ## Switching between different interaction function types
        potential_type_switch = {2:3,3:2,4:5,5:4}

        for p in range(self.n_model_param):
            p_pairs = self.model_param_interactions[p]
            for n in range(len(p_pairs)):
                if new_model_param_values[p] < 0.:
                    ## If parameter changes sign then the pairwise_type is flipped.
                    if self.pairwise_type[p_pairs[n]] == 1:
                        self.pairwise_type[p_pairs[n]] = 3
                    else:
                        self.potential_type[p_pairs[n]] = potential_type_switch[self.pairwise_type[p_pairs[n]]]
                else:
                    if self.pairwise_type[p_pairs[n]] == 1:
                        self.pairwise_type[p_pairs[n]] = 2
                ## Model parameters are always positive
                self.model_param_values[p] = abs(new_model_param_values[p])   

        ## Refresh everything that depends on model parameters
        self._set_nonbonded_interactions()

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
