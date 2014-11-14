""" SmogCalpha

Description:

    Generates topology and grofile needed to run Smog-style C-alpha Go-model. The
style of the input file formats emulates the SMOG server (see reference (1)).

Example Usage:
    See project_tools/examples

References:
(1) Noel, J.K.; Whitford, P.C.; Sanbonmatsu, K.Y.; Onuchic, J.N. SMOG@ctbp:
Simplified Deployment of Structure-Based Models in GROMACS. Nucleic Acids Res.
2010, 38, W657-61.
"""

import numpy as np
import subprocess as sb
import os
import time

import bonded_potentials as bond
import pairwise_potentials as pairwise
import pdb_parser


class SmogCalpha(object):
    """ This class creates a smog-like topology and grofile """
    ## To Do:
    ##  1. Delineate between 'model parameters' and 'interaction strengths':
    ##     where model parameters are non-redundant list.
    ##  2. Create list of interaction potentials.

    ## To Do:  Refactor argument passing to use *args or **kwargs:
    ##  def __init__(self,pdb,**kwargs):

    #def __init__(self,pdb,contacts=None,contact_epsilons=None,LJtype=None,contact_widths=None,noncontact_wall=None,
    #        epsilon_bar=None,contact_type=None,
    #        disulfides=None,model_code=None,contact_params=None,
    #        Tf_iteration=0,Mut_iteration=0,
    #        fitting_data=None,fitting_includes=[None],fitting_solver="Levenberg",fitting_allowswitch=False,
    #        dry_run=False):
    def __init__(self,**kwargs):

        self.path = os.getcwd()
        for key in kwargs.iterkeys():
            #print key.lower(),kwargs[key]
            if key in ["LJtype","Tf_iteration","Mut_iteration"]:
                setattr(self,key,kwargs[key])
            else:
                setattr(self,key.lower(),kwargs[key])

        if not os.path.exists(self.pdb):
            print "ERROR! The inputted pdb: %s does not exist" % pdb
            print " Exiting."
            raise SystemExit

        self.beadmodel = "CA"
        self.backbone_params = ["Kb","Ka","Kd"]
        self.backbone_param_vals = {"Kb":20000.,"Ka":400.,"Kd":1}
        self.citation = self.citation_info(self.model_code)
        self.error = 0
        self.initial_T_array = None
        self.name = self.pdb.split(".pdb")[0]
        self.subdir = self.name

        """
        self.model_code = model_code
        self.citation = self.citation_info(self.model_code)


        ## Contacts, their absolute strengths, attractive/repulsive character
        self.contact_type = contact_type
        self.contacts = contacts
        self.contact_params = contact_params
        self.contact_epsilons = contact_epsilons
        self.LJtype = LJtype
        self.contact_widths = contact_widths
        self.noncontact_wall = noncontact_wall
        
        ## Average contact strength
        self.epsilon_bar = epsilon_bar  

        self.disulfides = disulfides
        self.dry_run = dry_run
        self.error = 0
        self.initial_T_array = None


        self.Tf_iteration = Tf_iteration
        self.Mut_iteration = Mut_iteration
        
        ## Options for parameter fitting
        self.fitting_data = fitting_data
        self.fitting_includes = fitting_includes
        self.fitting_solver = fitting_solver
        self.fitting_allowswitch = fitting_allowswitch

        self.pdb = pdb
        self.name = pdb.split(".pdb")[0]
        self.subdir = self.name

        """
        ## Prepare the coordinates from the pdb file.

        self.cleanpdb_full, self.cleanpdb_full_noH, self.cleanpdb = pdb_parser.clean(self.pdb)
        self.dissect_native_pdb(self.cleanpdb)
        self.n_contacts = len(self.contacts)
        if self.LJtype == None:
            self.n_repcontacts = 0
        else:
            self.n_repcontacts = sum((self.LJtype == -1).astype(int))

        self.Qref = np.zeros((self.n_residues,self.n_residues))
        for pair in self.contacts:
            self.Qref[pair[0]-1,pair[1]-1] = 1 
        self._get_index_ndx()

        ## Check disulfide separation and remove from contacts list.
        self.check_disulfides() 

        ## Generate grofile, topology file, and necessary table files.
        self.check_contact_opts()
        self.generate_grofile()
        self.generate_topology()
        self._get_interaction_tables()
            
    def __repr__(self):
        ''' The string representation of all the model info.'''
        repstring = "[ Path ]\n"
        repstring += "%s\n" % self.path
        repstring += "[ PDB ]\n"
        repstring += "%s\n" % self.pdb
        repstring += "[ Subdir ]\n"
        repstring += "%s\n" % self.subdir
        repstring += "[ Tf_Iteration ]\n"
        repstring += "%s\n" % self.Tf_iteration
        repstring += "[ Mut_Iteration ]\n"
        repstring += "%s\n" % self.Mut_iteration
        repstring += "[ Model_Code ]\n" 
        repstring += "%s\n" % self.model_code
        repstring += "[ Bead_Model ]\n" 
        repstring += "%s\n" % self.beadmodel
        repstring += "[ Backbone_params ]\n" 
        repstring += "  %5s       %5s       %5s\n" % ("Kb","Ka","Kd")
        repstring += "[ Backbone_param_vals ]\n" 
        repstring += "%10.2f%10.2f%10.2f\n" % \
            (self.backbone_param_vals["Kb"],self.backbone_param_vals["Ka"],self.backbone_param_vals["Kd"])
        repstring += "[ Disulfides ]\n"
        if self.disulfides == None:
            repstring += "%s\n" % None
        else:
            temp = ''
            for x in self.disulfides:
                temp += " %d " % x 
            repstring += "%s\n" % temp
        repstring += "[ Epsilon_Bar ]\n"
        repstring += "%s\n" % str(self.epsilon_bar)
        repstring += "[ Contact_Params ]\n"
        repstring += "%s\n" % str(self.contact_params)
        repstring += "[ Contact_Type ]\n"
        repstring += "%s\n" % self.contact_type
        repstring += "[ Fitting_Data ]\n"
        repstring += "%s\n" % self.fitting_data
        repstring += "[ Fitting_Includes ]\n"
        for dir in self.fitting_includes:
            repstring += "%s" % str(dir)
        repstring += "\n"
        repstring += "[ Fitting_Solver ]\n"
        repstring += "%s\n" % self.fitting_solver
        repstring += "[ Fitting_AllowSwitch ]\n"
        repstring += "%s\n" % self.fitting_allowswitch
        repstring += "[ Reference ]\n" 
        repstring += "%s\n" % self.citation
        return repstring

    def append_log(self,string):
        now = time.localtime()
        now_string = "%s:%s:%s:%s:%s" % (now.tm_year,now.tm_mon,now.tm_mday,now.tm_hour,now.tm_min)
        logfile = open('%s/%s/%s.log' % (path,self.subdir,self.subdir),'a').write("%s %s\n" % (now_string,string))

    def citation_info(self,key):
        ''' Dictionary of references'''
        citations = {'HomGo':"Clementi, C.; Nymeyer, H.; Onuchic, J. N. \n"+  \
                    "Topological and Energetic Factors: What Determines the \n"+ \
                    "Structural Details of the Transition State Ensemble and \n"+ \
                    "'En-Route' Intermediates for Protein Folding? An Investigation\n"+ \
                    "for Small Globular Proteins. J. Mol. Biol. 2000, 298, 937-53.", 
                     'HetGo':"Matysiak, S; Clementi, C. Optimal Combination of \n"+ \
                    "Theory and Experiment for the Characterization of the Protein \n"+ \
                    "Folding Landscape of S6: How Far Can a Minimalist Model Go? \n"+ \
                    "J. Mol. Biol. 2004, 343, 235-48", 
                     'DMC':"Matyiak, S; Clementi, C. Minimalist Protein Model as \n"+  \
                    "a Diagnostic Tool for Misfolding and Aggregation. J. Mol. Biol. \n"+ \
                    "2006, 363, 297-308.",None:"None"}
        return citations[key]

    def check_contact_opts(self):
        """ """
        if self.contact_epsilons == None:
            print "  No contact epsilons set. Setting contacts to homogeneous model. 1"
            self.contact_epsilons = np.ones(self.n_contacts,float)
        if self.contact_type == "LJ1210":
            if self.LJtype == None:
                print "  No LJtype given. Setting contacts to attractive. 1"
                self.LJtype = np.ones(self.n_contacts,float)
        elif self.contact_type == "Gaussian":
            if self.contact_widths == None:
                print "  No contact widths set. Setting contact widths to attractive. 0.5 Angstrom"
                self.contact_widths = 0.05*np.ones(self.n_contacts,float)
            if self.noncontact_wall == None:
                print "  No noncontact wall set. Setting noncontact wall to 0.4 nm"
                noncontact = 0.4
                self.noncontact_wall = noncontact*np.ones(self.n_contacts,float)

        self.interaction_potentials = []

    def check_disulfides(self):
        ''' Check that specified disulfides are between cysteine and that 
            the corresonding pairs are within 0.8 nm.

        To Do: 
            - Add disulfides to the list of atoms to be bonded.
        '''
        coords = self.atom_coords
        residues = self.atom_residues
        if self.disulfides != None:
            print "  Checking disulfides are reasonable."
            for i in range(len(self.disulfides[::2])):
                cys1 = self.disulfides[2*i]
                cys2 = self.disulfides[2*i + 1]
                if (residues[cys1-1] != "CYS") or (residues[cys2-1] != "CYS"):
                    print "WARNING! Specifying disulfide between two residues that aren't CYS cysteine! "
                    #print "Exiting"
                    #raise SystemExit
                #print "##  Checking disulfide residue identities: ", residues[cys1-1], residues[cys2-1]    ## DEBUGGING
                dist = coords[cys1-1] - coords[cys2-1]
                separation = np.linalg.norm(dist)
                #print "##  Separation: ", separation ## DEBUGGING
                if separation > 0.8:
                    print "WARNING! Specifying disulfide of separation greater than 0.8 nm."
                    #print "Specifying disulfide of separation greater than 0.8 nm."
                    #print "Exiting"
                    #raise SystemExit
                else:
                    print "   ",residues[cys1-1]+str(cys1), residues[cys2-1]+str(cys2), \
                          " are separated by: %.4f nm. Good." % separation
                    ## Remove disulfide pair from self.contacts if it is there.
                    new_conts = []
                    for pair in self.contacts:
                        if (pair[0] == cys1) and (pair[1] == cys2):
                            continue
                        else:
                            new_conts.append(pair)
                    self.contacts = np.array(new_conts)

                if self.Qref[cys1-1][cys2-1] == 1:
                    #print "    Subtracting 1 from n_contacts for ", cys1,cys2, " disulfide"
                    self.n_contacts -= 1
        else:
            print "  No disulfides to check."
        self.contacts_ndx = "[ contacts ]\n"
        for i in range(self.n_contacts):
            self.contacts_ndx += "%4d %4d\n" % (self.contacts[i][0],self.contacts[i][1])


    def dissect_native_pdb(self,pdb):
        ''' Extract info from the Native.pdb for making index and top file'''
        indices = []
        atoms = []
        residues = []
        coords = []
        pdblines = pdb.split("\n")
        res_indx = 0
        for line in pdblines:
            if line.startswith("END"):
                break
            else:
                indices.append(int(line[6:13]))
                atoms.append(line[11:16].strip())
                coords.append([float(line[31:39]),float(line[39:47]),float(line[47:55])]) 
                if (int(line[23:26]) == (res_indx + 1)):
                    res_indx += 1 
                    residues.append(line[17:20])

        ## Coordinates in pdb files are Angstroms. Convert to nanometers.
        coords = np.array(coords)/10.
        self.n_residues = len(residues)
        self.n_atoms = len(atoms)
        self.atom_indices = indices
        self.atom_types = atoms
        self.atom_residues = residues
        self.atom_coords = coords

        ## Set bonded force field terms quantities.
        self.bonded_indices = [[indices[i],indices[i+1]] for i in range(self.n_atoms-1)]
        self.angled_indices = [[indices[i],indices[i+1],indices[i+2]] for i in range(self.n_atoms-2)]
        self.dihedral_indices = [[indices[i],indices[i+1],indices[i+2],indices[i+3]] for i in range(self.n_atoms-3)]

        self.bond_min = [ bond.distance(coords,i_idx-1,j_idx-1) for i_idx,j_idx in self.bonded_indices ]
        self.angle_min = [ bond.angle(coords,i_idx-1,j_idx-1,k_idx-1) for i_idx,j_idx,k_idx in self.angled_indices ]
        self.dihedral_min = [ bond.dihedral(coords,i_idx-1,j_idx-1,k_idx-1,l_idx-1) for i_idx,j_idx,k_idx,l_idx in self.dihedral_indices ]

        if not hasattr(self,"bond_strengths"):
            self.bond_strengths = [ self.backbone_param_vals["Kb"] for i in range(len(self.bond_min)) ]
        if not hasattr(self,"angle_strengths"):
            self.angle_strengths = [ self.backbone_param_vals["Ka"] for i in range(len(self.angle_min)) ]
        if not hasattr(self,"dihedral_strengths"):
            self.dihedral_strengths = [ self.backbone_param_vals["Kd"] for i in range(len(self.dihedral_min)) ]


    def calculate_contact_potential(self,rij):
        """ Contact potential for all contacts in the trajectory """

        ## To Do:
        ## 1. Loop over 

        ## Epsilons, deltas, and sigmas for all contacts
        conts = np.arange(self.n_contacts)
        epsilons = self.contact_epsilons[conts]
        sigmas = self.contact_sigmas[conts]
        if self.contact_type == "LJ1210":
            allindx = np.arange(self.n_contacts)
            repindx = allindx[self.LJtype == -1]
            attindx = allindx[self.LJtype == 1]

            ## Only count values of potential energy function where interaction is
            ## attractive. Assign attractive interaction en masse then loop over
            ## repulsive interactions because each one requires different function
            ## evaluation.
            Vij = np.zeros(rij.shape,float)
            x = sigmas/rij
            #x[x > 1.09] = 1.09  # <-- 1.09 is where LJ12-10 crosses zero for attractive
            Vij[:,attindx] = epsilons[attindx]*(5.*(x[:,attindx]**12) - 6.*(x[:,attindx]**10))     
            for i in range(self.n_repcontacts):
                indx = repindx[i]
                eps = self.contact_epsilons[indx]
                sigma = self.contact_sigmas[indx]
                x_indx = x[:,indx]
                V_indx = self._get_repLJ1210_potential(x_indx,eps,sigma)
                Vij[:,indx] = V_indx

        elif self.contact_type == "Gaussian":
            ## NEEDS TO BE TESTED
            noncontact_sigmas = self.noncontact_wall[conts]
            contact_sigmas = np.ones(rij.shape)*sigmas
            widths = self.contact_widths[conts]
            contact_widths = np.ones(rij.shape)*widths
            
            x = noncontact_sigmas/rij
            Vij = epsilons*((1. + (1./epsilons)*(x**12))*(1. - np.exp(-((rij - contact_sigmas)**2)/(2.*(contact_widths**2)))) - 1.)

        return Vij


    def _get_interaction_tables(self):
        ''' Returns the table for interaction type 'i'. The values of the
            potential is set to 0.0 when r is too close to zero to avoid
            blowup. 
        '''

        if self.contact_type == "LJ1210":
    
            self.table = self._get_LJ1210_table()

            self.rep_tables = []
            self.rep_tablenames = []
            if self.n_repcontacts != 0:
                for i in range(self.n_repcontacts):
                    pair = self.contacts[self.LJtype == -1][i]
                    epsilon = self.contact_epsilons[self.LJtype == -1][i]
                    sigma = self.contact_sigmas[self.LJtype == -1][i]
                    table = self._get_repLJ1210_table(epsilon,sigma)
                    table_name = "table_b%d.xvg" % (i+1)

                    self.rep_tables.append(table)
                    self.rep_tablenames.append(table_name)

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

    def _get_repLJ1210_table(self,eps,sigma):
        ''' A repulsive LJ12-10 interaction potential '''
        r = np.arange(0.0,100.0,0.002)
        r[0] = 1
        x = sigma/r
        table = np.zeros((len(r),3),float)
        table[:,0] = r
        table[x > 1,1] = eps*(5*(x[x > 1]**12) - 6*(x[x > 1]**10) + 2)
        table[x <= 1,1] = -eps*(5*(x[x <= 1]**12) - 6*(x[x <= 1]**10))
        table[x > 1, 2] = 60*(eps/sigma)*(x[x > 1]**13 - x[x > 1]**11)
        table[x <= 1,2] = -60*(eps/sigma)*(x[x <= 1]**13 - x[x <= 1]**11)
        table[:5, 1:] = 0.
        table[0,0] = 0

        return table

    def _get_repLJ1210_potential(self,x,eps,sigma):
        V = np.zeros(x.shape,float)
        V[x > 1] = eps*(5*(x[x > 1]**12) - 6*(x[x > 1]**10) + 2)
        V[x <= 1] = -eps*(5*(x[x <= 1]**12) - 6*(x[x <= 1]**10))
        return V

    def _get_index_ndx(self):
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
            i_idx = self.bonded_indices[j][0]
            j_idx = self.bonded_indices[j][1]
            dist = self.bond_min[j]
            kb = self.bond_strengths[j]
            bonds_string += "%6d %6d%2d%18.9e%18.9e\n" %  \
                          (i_idx,j_idx,1,dist,kb)
        return bonds_string

    def _get_tabled_string(self):
        tabled_string = ""
        ## Add special nonbonded table interactions. 
        if (self.contact_type == "LJ1210"):
            if self.n_repcontacts != 0:
                tabled_string += "; tabled interactions contacts below\n"
                for i in range(self.n_repcontacts):
                    ## add extra tabled string
                    pair = self.contacts[self.LJtype == -1][i]
                    tabled_string += "%6d %6d%2d%18d%18.9e\n" %  \
                              (pair[0],pair[1],9,i+1,1.0)

        return tabled_string

    def _get_angles_string(self):
        ''' Generate the [ angles ] string.'''
        angles_string = " [ angles ]\n"
        angles_string += " ; ai  aj  ak  func  th0(deg)   Ka\n"
        for j in range(len(self.angle_min)):
            i_idx = self.angled_indices[j][0]
            j_idx = self.angled_indices[j][1]
            k_idx = self.angled_indices[j][2]
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
        """ Get the [ pairs ] string"""

        if self.epsilon_bar != None:
            print "  Scaling attractive contacts such that epsilon_bar =", self.epsilon_bar
            ## TO DO:
            #stability = sum(self.contact_epsilons[self.LJtype == 1]) - sum(
            #avg_eps = np.mean(self.contact_epsilons[attractive])
            #self.contact_epsilons[attractive] = /avg_eps

        self.contact_sigmas = np.zeros(len(self.contacts),float)
        pairs_string = " [ pairs ]\n"
        pairs_string += " ; i j type and weights\n"
        beadbead_string = ""
        for i in range(len(self.contacts)):
            
            res_a = self.contacts[i][0]
            res_b = self.contacts[i][1]
            x_a = self.atom_coords[res_a-1]
            x_b = self.atom_coords[res_b-1]
            sig_ab = np.linalg.norm(x_a - x_b)
            self.contact_sigmas[i] = sig_ab
            eps_ab = self.contact_epsilons[i] 

            ## Generalize for writing parameters for any pair potential
            if self.contact_type in [None, "LJ1210"]:
                if self.LJtype[i] == 1:
                    c12 = eps_ab*5.0*(sig_ab**12)
                    c10 = eps_ab*6.0*(sig_ab**10)

                    pairs_string += "%6d %6d%2d%18.9e%18.9e\n" % \
                                    (res_a,res_b,1,c10,c12)

                resid_a = self.atom_residues[res_a-1]
                resid_b = self.atom_residues[res_b-1]

                beadbead_string += "%5d%5d%8s%8s%5d%18.9e%18.9e%18d\n" % \
                                (res_a,res_b,resid_a,resid_b,i+1,sig_ab,eps_ab,self.LJtype[i])
            elif self.contact_type == "Gaussian":
                width_ab = self.contact_widths[i] 
                ncwall_ab = self.noncontact_wall[i]**12
                pairs_string += "%6d %6d%2d%18.9e%18.9e%18.9e%18.9e\n" % \
                                (res_a,res_b,6,eps_ab,sig_ab,width_ab,ncwall_ab)

                resid_a = self.atom_residues[res_a-1]
                resid_b = self.atom_residues[res_b-1]

                beadbead_string += "%5d%5d%8s%8s%5d%18.9e%18.9e%18d%18.9e%18.9e\n" % \
                    (res_a,res_b,resid_a,resid_b,i+1,sig_ab,eps_ab,1,width_ab,ncwall_ab)
            else:
                print "ERROR! unrecognized contact option."
                print "Exiting"
                raise SystemExit

        noncontact = 0.4        ## Noncontact radius 4 Angstroms
        if self.disulfides != None:
            pairs_string, beadbead_string = self._add_disulfides(pairs_string,beadbead_string,noncontact)
        self.beadbead = beadbead_string 
        return pairs_string

    def _add_disulfides(self,pairs_string,beadbead_string,noncontact):
        ## To Do: Add disulfides to list of bonded interactions.
        for i in range(len(self.disulfides)/2):
            cys_a = self.disulfides[2*i]
            cys_b = self.disulfides[2*i + 1]
            x_a = self.atom_coords[cys_a-1]
            x_b = self.atom_coords[cys_b-1]
            sig_ab = np.linalg.norm(x_a - x_b)
            eps_ab = 50.
            width_ab = 0.05
            noncontact = 0.4**12

            if self.contact_type in [None, "LJ1210"]:
                c12 = eps_ab*5.0*(sig_ab**12)
                c10 = eps_ab*6.0*(sig_ab**10)

                print " Linking disulfide: ",cys_a,cys_b, " with eps = ",eps_ab
                pairs_string += "%6d %6d%2d%18.9e%18.9e\n" % \
                                (cys_a,cys_b,1,c10,c12)

                resid_a = self.atom_residues[cys_a-1]
                resid_b = self.atom_residues[cys_b-1]

                beadbead_string += "%5d%5d%8s%8s%5s%18.9e%18.9e%18d\n" % \
                                (cys_a,cys_b,resid_a,resid_b,"ss",sig_ab,eps_ab,1)
            elif self.contact_type == "Gaussian":
                print " Linking disulfide: ",cys_a,cys_b, " with eps = ",eps_ab
                pairs_string += "%6d %6d%2d%18.9e%18.9e%18.9e%18.9e\n" % \
                                (res_a,res_b,6,eps_ab,sig_ab,width_ab,noncontact)

                resid_a = self.atom_residues[res_a-1]
                resid_b = self.atom_residues[res_b-1]

                beadbead_string += "%5d%5d%8s%8s%5d%18.9e%18.9e%18d%18.9e%18.9e\n" % \
                    (res_a,res_b,resid_a,resid_b,"ss",sig_ab,eps_ab,1,width_ab,noncontact)

        return pairs_string, beadbead_string


    def _get_exclusions_string(self):
        """ Get the [ exclusions ] string"""
        exclusions_string = " [ exclusions ]\n"
        exclusions_string += " ; ai aj \n"
        for i in range(len(self.contacts)):
            res_a = self.contacts[i][0]
            res_b = self.contacts[i][1]
            exclusions_string += "%6d %6d\n" % (res_a,res_b)
        if self.disulfides != None:
            for i in range(len(self.disulfides)/2):
                cys_a = self.disulfides[2*i]
                cys_b = self.disulfides[2*i + 1]
                exclusions_string += "%6d %6d\n" % (cys_a,cys_b)

        return exclusions_string

    def generate_grofile(self):

        gro_string = " Structure-based Gro file\n"
        gro_string += "%12d\n" % len(self.atom_types)
        for i in range(len(self.atom_types)):
            gro_string += "%5d%5s%4s%6d%8.3f%8.3f%8.3f\n" % \
                (self.atom_indices[i],self.atom_residues[i],self.atom_types[i],self.atom_indices[i],
                self.atom_coords[i][0],self.atom_coords[i][1],self.atom_coords[i][2])
        gro_string += "   %-25.16f%-25.16f%-25.16f" % (5.0,5.0,5.0)

        self.grofile = gro_string

    def generate_topology(self):
        """ Return a structure-based topology file. SMOG-style."""

        top_string =  " ; Structure-based  topology file for Gromacs:\n"
        top_string += " [ defaults ]\n"
        top_string += " ;nbfunc comb-rule gen-pairs\n"
        top_string += "      1           1 no\n\n"
        top_string += " [ atomtypes ]\n"
        top_string += " ;name  mass     charge   ptype c10       c12\n"
        top_string += " CA     1.000    0.000 A    0.000   0.167772160E-04\n\n"
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

    model = SmogCalpha(pdb,contacts=contacts,contact_type="Gaussian")

    open("%s.top" % name, "w").write(model.topology)
    open("%s.gro" % name, "w").write(model.grofile)

