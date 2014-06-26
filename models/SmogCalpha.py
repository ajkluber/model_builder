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

class SmogCalpha(object):
    """ This class creates a smog-like topology and grofile """

    def __init__(self,pdb,contacts=None,contact_epsilons=None,contact_deltas=None,
            epsilon_bar=None,disulfides=None,modelcode=None,contact_energies=None,
            Tf_iteration=0,Mut_iteration=0,dryrun=False):

        self.beadmodel = "CA"
        self.backbone_params = ["Kb","Ka","Kd"]
        self.backbone_param_vals = {"Kb":20000.,"Ka":400.,"Kd":1}

        self.modelcode = modelcode
        self.citation = self.citation_info(self.modelcode)
        self.contact_energies = contact_energies
        self.contact_epsilons = contact_epsilons
        self.contact_deltas = contact_deltas
        self.epsilon_bar = epsilon_bar
        self.disulfides = disulfides
        self.dryrun = dryrun
        self.error = 0
        self.initial_T_array = None

        self.Tf_iteration = Tf_iteration
        self.Mut_iteration = Mut_iteration

        self.path = os.getcwd()

        if not os.path.exists(pdb):
            print "ERROR! The inputted pdb: ",pdb," does not exist"
            print " Exiting."
            raise SystemExit
        else:
            self.pdb = pdb
            self.name = pdb.split(".pdb")[0]
            self.subdir = self.name
            if not os.path.exists(self.subdir+"/Qref_shadow"):
                os.makedirs(self.subdir+"/Qref_shadow")

            self.clean_pdb()
            self.dissect_native_pdb()
            self.get_index_ndx()
            if contacts != None:
                self.skip_shadow_contacts(contacts)
            else:
                self.shadow_contacts()
            self.check_disulfides() 
            self.generate_grofile()
            self.generate_topology()
            self.get_interaction_table()
            
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
        repstring += "%s\n" % self.modelcode
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
        repstring += "[ Contact_Energies ]\n"
        repstring += "%s\n" % str(self.contact_energies)
        repstring += "[ Reference ]\n" 
        repstring += "%s\n" % self.citation
        return repstring

    def append_log(self,string):
        now = time.localtime()
        now_string = "%s:%s:%s:%s:%s" % (now.tm_year,now.tm_mon,now.tm_mday,now.tm_hour,now.tm_min)
        logfile = open(self.path+'/'+self.subdir+'/'+self.subdir+'.log','a').write(now_string+' '+string+'\n')

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

    def check_disulfides(self):
        ''' Check that specified disulfides are between cysteine and that 
            the corresonding pairs are within 0.8 nm.'''
        coords = self.coords
        residues = self.residues
        if self.disulfides != None:
            print "  Checking disulfides are reasonable."
            for i in range(len(self.disulfides[::2])):
                cys1 = self.disulfides[2*i]
                cys2 = self.disulfides[2*i + 1]
                if (residues[cys1-1] != "CYS") or (residues[cys2-1] != "CYS"):
                    print "ERROR! Specifying disulfide between two residues that aren't CYS cysteine! "
                    print "Exiting"
                    raise SystemExit
                #print "##  Checking disulfide residue identities: ", residues[cys1-1], residues[cys2-1]    ## DEBUGGING
                dist = coords[cys1-1] - coords[cys2-1]
                separation = np.linalg.norm(dist)
                #print "##  Separation: ", separation ## DEBUGGING
                if separation > 0.8:
                    print "ERROR!"
                    print "Specifying disulfide of separation greater than 0.8 nm."
                    print "Exiting"
                    raise SystemExit
                else:
                    print "   ",residues[cys1-1]+str(cys1), residues[cys2-1]+str(cys2), \
                          " are separated by: %.4f nm. Good." % separation
                    ## Remove disulfide pair from self.contacts if it is there.
                    new_conts = []
                    for pair in self.contacts:
                        if all(pair == [cys1,cys2]):
                            pass
                        else:
                            new_conts.append(pair)
                    self.contacts = np.array(new_conts)
                if self.Qref[cys1-1][cys2-1] == 1:
                    #print "    Subtracting 1 from n_contacts for ", cys1,cys2, " disulfide"
                    self.n_contacts -= 1
        else:
            print "  No disulfides to check."

    def clean_pdb(self):
        """ Grab only the lines of the pdb that we want. 
    
        Description:
            
            Returns an all-atom pdb and a C-alpha only pdb.

        Find the full PDB format specification at:
        http://www.wwpdb.org/documentation/format33/sect9.html#ATOM

        PDB fixed-width column format is given by:
        ATOM     44  C   ALA A  11      12.266  21.667  20.517  1.00 28.80           C  
        """
        first_full = 0
        atomid_full = 1
        cleanpdb_full = ''
        cleanpdb_full_noH = ''
        first_ca = 0
        atomid_ca = 1
        cleanpdb_ca = ''
        for line in open(self.pdb,'r'):
            line = line.rstrip("\n")
            if line[:3] in ['TER','END']:
                break
            else:
                ## Keep only ATOM lines.
                if line[:4] == 'ATOM':
                    if line[13:16].strip() == "CA":
                        if first_ca == 0:
                            if line[16] in ["A"," "]:
                                newline_ca = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                        (atomid_ca,line[12:16],line[17:20],1,line[26:55])
                                atomid_ca += 1
                                first_ca = 1
                                first_index_ca = int(line[22:26]) - 1
                                cleanpdb_ca += newline_ca
                        else:
                            if line[16] in ["A"," "]:
                                newline_ca = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                        (atomid_ca,line[12:16],line[17:20],int(line[22:26])-first_index_ca,line[26:55])
                                atomid_ca += 1
                                cleanpdb_ca += newline_ca

                    if first_full == 0:
                        if (line[16] in ["A"," "]) and (line[13] not in ["E","D"]):
                            newline_full = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid_full,line[12:16],line[17:20],1,line[26:55])
                            atomid_full += 1
                            first_full = 1
                            first_index_full = int(line[22:26]) - 1
                            cleanpdb_full += newline_full
                            ## strip Hydrogens
                            if not line[12:16].strip().startswith("H"):
                                cleanpdb_full_noH += newline_full
                    else:
                        if (line[16] in ["A"," "]) and line[13] not in ["E","D"]:
                            newline_full = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid_full,line[12:16],line[17:20],int(line[22:26])-first_index_full,line[26:55])
                            atomid_full += 1
                            cleanpdb_full += newline_full
                            if not line[12:16].strip().startswith("H"):
                                cleanpdb_full_noH += newline_full
         
        cleanpdb_full += 'END\n'
        cleanpdb_full_noH += 'END\n'
        cleanpdb_ca += 'END\n'
        self.cleanpdb = cleanpdb_ca
        self.cleanpdb_full = cleanpdb_full
        self.cleanpdb_full_noH = cleanpdb_full_noH
        open(self.subdir+"/Native.pdb","w").write(self.cleanpdb)
        open(self.subdir+"/Qref_shadow/clean.pdb","w").write(self.cleanpdb_full)
        open(self.subdir+"/clean.pdb","w").write(self.cleanpdb_full)
        open(self.subdir+"/clean_noH.pdb","w").write(self.cleanpdb_full_noH)
        open(self.subdir+"/"+self.subdir+".pdb","w").write(self.cleanpdb_full_noH)

    def dissect_native_pdb(self):
        ''' Extract info from the Native.pdb for making index and top file'''
        indices = []
        atoms = []
        residues = []
        coords = []

        for line in open(self.subdir+"/Native.pdb","r"):
            if line.startswith("END"):
                break
            else:
                indices.append(int(line[6:13]))
                atoms.append(line[11:16].strip())
                residues.append(line[17:20])
                coords.append([float(line[31:39]),float(line[39:47]),float(line[47:55])]) 

        ## Coordinates in pdb files are Angstroms. Convert to nanometers.
        coords = np.array(coords)/10.
        self.indices = indices
        self.atoms = atoms
        self.residues = residues
        self.coords = coords
        self.n_residues = len(residues)

    def dihedral(self,coords,i_idx,j_idx,k_idx,l_idx):
        ''' Compute the dihedral between planes. '''
        v21 = coords[j_idx] - coords[i_idx]
        v31 = coords[k_idx] - coords[i_idx]
        v32 = coords[k_idx] - coords[j_idx]
        v42 = coords[l_idx] - coords[j_idx]
        v21xv31 = np.cross(v21,v31)
        v21xv31 /= np.linalg.norm(v21xv31)
        v32xv42 = np.cross(v32,v42)
        v32xv42 /= np.linalg.norm(v32xv42)
        if np.dot(v21,v32xv42) < 0.:
            sign = -1.
        else:
            sign = 1.
        phi = 180. + sign*(180./np.pi)*np.arccos(np.dot(v21xv31,v32xv42))
        return phi

    def get_interaction_table(self):
        ''' Returns the table for interaction type 'i'. The values of the
            potential is set to 0.0 when r is too close to zero to avoid
            blowup. 
        '''
        r = np.arange(0.0,4.0,0.002)
        self.table = np.zeros((len(r),7),float)
        self.table[:,0] = r
        self.table[20:,1:7] = self.get_LJ1210_table(r[20:])

    def get_LJ1210_table(self,r):
        ''' LJ12-10 interaction potential ''' 
        table = np.zeros((len(r),6),float)
        table[:,0] = 1./r
        table[:,1] = 1./(r**2)
        table[:,2] = -1./(r**10)
        table[:,3] = -10./(r**11)
        table[:,4] = 1./(r**12)
        table[:,5] = 12./(r**13)
        return table

    def get_residue_mass(self):
        '''Masses in atomic mass units.'''
        residue_mass = {'ALA':   89.0935, 'ARG':  174.2017, 'ASN':  132.1184,
                        'ASP':  133.1032, 'CYS':  121.1590, 'GLN':  146.1451,
                        'GLU':  147.1299, 'GLY':   75.0669, 'HIS':  155.1552,
                        'ILE':  131.1736, 'LEU':  131.1736, 'LYS':  146.1882,
                        'MET':  149.2124, 'PHE':  165.1900, 'PRO':  115.1310,
                        'SER':  105.0930, 'THR':  119.1197, 'TRP':  204.2262,
                        'TYR':  181.1894, 'VAL':  117.1469, 'SOL':   18.0150}
        return residue_mass

    def get_residue_one_letter_code(self):
        '''Converting from three letter code to one letter FASTA code.'''
        residue_code = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N',
                        'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
                        'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
                        'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                        'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                        'SER': 'S', 'THR': 'T', 'TRP': 'W',
                        'TYR': 'Y', 'VAL': 'V'}
        return residue_code

    def get_index_ndx(self):
        ''' Generates index file for gromacs analysis utilities. '''
        ca_string = ''
        i = 1
        
        for indx in self.indices: 
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

    def get_atoms_string(self):
        ''' Generate the [ atoms ] string.'''
        atoms_string = ""
        for j in range(len(self.indices)):
            atomnum = self.indices[j]
            resnum = self.indices[j]
            resid = self.residues[j]
            atoms_string += " %5d%4s%8d%5s%4s%8d%8.3f%8.3f\n" % \
                        (atomnum,"CA",resnum,resid,"CA",atomnum,0.0,1.0)
        return atoms_string

    def get_bonds_string(self):
        ''' Generate the [ bonds ] string.'''
        kb = self.backbone_param_vals["Kb"]
        bonds_string = ""
        for j in range(len(self.indices)-1):
            i_idx = self.indices[j]-1
            j_idx = self.indices[j+1]-1
            dist = np.linalg.norm(self.coords[i_idx] - self.coords[j_idx])
            bonds_string += "%6d %6d%2d%18.9e%18.9e\n" %  \
                          (i_idx+1,j_idx+1,1,dist,kb)
        return bonds_string

    def get_angles_string(self):
        ''' Generate the [ angles ] string.'''
        ka = self.backbone_param_vals["Ka"]
        angles_string = ""
        for j in range(len(self.indices)-2):
            i_idx = self.indices[j]-1
            j_idx = self.indices[j+1]-1
            k_idx = self.indices[j+2]-1
            xkj = self.coords[k_idx] - self.coords[j_idx]
            xkj /= np.linalg.norm(xkj)
            xij = self.coords[i_idx] - self.coords[j_idx]
            xij /= np.linalg.norm(xij)
            theta = (180./np.pi)*np.arccos(np.dot(xkj, xij))
            angles_string += "%6d %6d %6d%2d%18.9e%18.9e\n" %  \
                          (i_idx+1,j_idx+1,k_idx+1,1,theta,ka)
        return angles_string

    def get_dihedrals_string(self):
        ''' Generate the [ dihedrals ] string.'''
        kd = self.backbone_param_vals["Kd"]
        dihedrals_string = ""
        dihedrals_ndx = '[ dihedrals ]\n'
        for j in range(len(self.indices)-3):
            i_idx = self.indices[j]-1
            j_idx = self.indices[j+1]-1
            k_idx = self.indices[j+2]-1
            l_idx = self.indices[j+3]-1
            phi = self.dihedral(self.coords,i_idx,j_idx,k_idx,l_idx)
            dihedrals_string += "%6d %6d %6d %6d%2d%18.9e%18.9e%2d\n" %  \
                          (i_idx+1,j_idx+1,k_idx+1,l_idx+1,1,phi,kd,1)
            dihedrals_string += "%6d %6d %6d %6d%2d%18.9e%18.9e%2d\n" %  \
                          (i_idx+1,j_idx+1,k_idx+1,l_idx+1,1,3.*phi,kd/2.,3)
            dihedrals_ndx += '%4d %4d %4d %4d\n' % \
                                (i_idx+1,j_idx+1,k_idx+1,l_idx+1)
        self.dihedrals_ndx = dihedrals_ndx
        return dihedrals_string

    def get_pairs_string(self):
        """ Get the [ pairs ] string"""
        if self.contact_epsilons == None:
            print "  No contact epsilons set. Setting contacts to homogeneous model. 1"
            self.contact_epsilons = np.ones(len(self.contacts),float)
        if self.contact_deltas == None:
            print "  No contact deltas set. Setting contacts to attractive. 1"
            self.contact_deltas = np.ones(len(self.contacts),float)
        if self.epsilon_bar != None:
            print "  Scaling attractive contacts such that epsilon_bar =", self.epsilon_bar
            avg_eps = np.sum(self.contact_epsilons*self.contact_deltas)
            attractive = (self.contact_deltas == 1.)
            self.contact_epsilons[attractive] = self.contact_epsilons[attractive]*self.epsilon_bar/avg_eps

        pairs_string = ""
        beadbead_string = ""
        for i in range(len(self.contacts)):
            res_a = self.contacts[i][0]
            res_b = self.contacts[i][1]
            x_a = self.coords[res_a-1]
            x_b = self.coords[res_b-1]
            sig_ab = np.linalg.norm(x_a - x_b)
            eps_ab = self.contact_epsilons[i] 
            delta_ab = self.contact_deltas[i] 

            c12 = eps_ab*5.0*(sig_ab**12)
            c10 = delta_ab*6.0*(sig_ab**10)

            pairs_string += "%6d %6d%2d%18.9e%18.9e\n" % \
                            (res_a,res_b,1,c10,c12)

            resid_a = self.residues[res_a-1]
            resid_b = self.residues[res_b-1]

            beadbead_string += "%5d%5d%8s%8s%5d%18.9e%18.9e%18.9e\n" % \
                            (res_a,res_b,resid_a,resid_b,i+1,sig_ab,eps_ab,delta_ab)

        if self.disulfides != None:
            for i in range(len(self.disulfides)/2):
                cys_a = self.disulfides[2*i]
                cys_b = self.disulfides[2*i + 1]
                x_a = self.coords[cys_a-1]
                x_b = self.coords[cys_b-1]
                sig_ab = np.linalg.norm(x_a - x_b)
                eps_ab = 100.
                delta_ab = 1.

                c12 = eps_ab*5.0*(sig_ab**12)
                c10 = delta_ab*6.0*(sig_ab**10)

                print " Linking disulfide: ",cys_a,cys_b, " with eps = ",eps_ab
                pairs_string += "%6d %6d%2d%18.9e%18.9e\n" % \
                                (cys_a,cys_b,1,c10,c12)

                resid_a = self.residues[cys_a-1]
                resid_b = self.residues[cys_b-1]

                beadbead_string += "%5d%5d%8s%8s%5s%18.9e%18.9e%18.9e\n" % \
                                (cys_a,cys_b,resid_a,resid_b,"ss",sig_ab,eps_ab,delta_ab)

        self.beadbead = beadbead_string 
        return pairs_string

    def get_exclusions_string(self):
        """ Get the [ exclusions ] string"""

        exclusions_string = ""
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
        gro_string += "%12d\n" % len(self.atoms)
        for i in range(len(self.atoms)):
            gro_string += "%5d%5s%4s%6d%8.3f%8.3f%8.3f\n" % \
                (self.indices[i],self.residues[i],self.atoms[i],self.indices[i],
                self.coords[i][0],self.coords[i][1],self.coords[i][2])
        gro_string += "   %-25.16f%-25.16f%-25.16f" % (5.0,5.0,5.0)

        self.grofile = gro_string

    def generate_topology(self):
        """ Return a structure-based topology file. SMOG-style."""

        top_string =  " ; Structure-based  topology file for Gromacs:\n"
        top_string += " [ defaults ]\n"
        top_string += " ;nbfunc comb-rule gen-pairs\n"
        top_string += "      1           1 no\n"
        top_string += "\n"

        top_string += " [ atomtypes ]\n"
        top_string += " ;name  mass     charge   ptype c10       c12\n"
        top_string += " CA     1.000    0.000 A    0.000   0.167772160E-04\n"
        top_string += "\n"

        top_string += " [ moleculetype ]\n"
        top_string += " ;name   nrexcl\n"
        top_string += " Macromolecule           3\n"
        top_string += "\n"

        top_string += " [ atoms ]\n"
        top_string += " ;nr  type  resnr residue atom  cgnr charge  mass\n"
        atoms_string = self.get_atoms_string()
        top_string += atoms_string
        top_string += "\n"

        top_string += " [ bonds ]\n"
        top_string += " ; ai aj func r0(nm) Kb\n"
        bonds_string = self.get_bonds_string()
        top_string += bonds_string
        top_string += "\n"

        top_string += " [ angles ]\n"
        top_string += " ; ai  aj  ak  func  th0(deg)   Ka\n"
        angles_string = self.get_angles_string()
        top_string += angles_string
        top_string += "\n"

        top_string += " [ dihedrals ]\n"
        top_string += " ; ai  aj  ak al  func  phi0(deg)   Kd mult\n"
        dihedrals_string = self.get_dihedrals_string()
        top_string += dihedrals_string
        top_string += "\n"

        top_string += " [ pairs ]\n"
        top_string += " ; i j type and weights\n"
        pairs_string = self.get_pairs_string()
        top_string += pairs_string
        top_string += "\n"

        top_string += " [ exclusions ]\n"
        top_string += " ; ai aj \n"
        exclusions_string = self.get_exclusions_string()
        top_string += exclusions_string
        top_string += "\n"

        top_string += " [ system ]\n"
        top_string += " ; name\n"
        top_string += " Macromolecule\n"
        top_string += "\n"

        top_string += " [ molecules ]\n"
        top_string += " ; name molec \n"
        top_string += " Macromolecule 1\n"
        top_string += "\n"

        self.topology = top_string

    def skip_shadow_contacts(self,contacts):
        """ When contacts are an input. """
        cwd = os.getcwd()
        self.contacts = contacts
        self.n_contacts = len(self.contacts)
        self.Qref = np.zeros((self.n_residues,self.n_residues))
        for pair in contacts:
            self.Qref[pair[0]-1,pair[1]-1] = 1 

        np.savetxt(cwd+"/"+self.subdir+"/Qref_shadow/Qref_cryst.dat",self.Qref,delimiter=" ",fmt="%1d")
        np.savetxt(cwd+"/"+self.subdir+"/Qref_cryst.dat",self.Qref,delimiter=" ",fmt="%1d")
        np.savetxt(cwd+"/"+self.subdir+"/contacts.dat",self.contacts,delimiter=" ",fmt="%4d")
        self.contacts_ndx = "[ contacts ]\n"
        for i in range(self.n_contacts):
            self.contacts_ndx += "%4d %4d\n" % (self.contacts[i][0],self.contacts[i][1])

    def shadow_contacts(self):
        ''' Call SMOG Shadow jar code to determine the shadow contacts. If 
            the reference matrix Qref_cryst.dat doesn't exist then create 
            and dive into a subdirectory to run shadow map. Then save 
            Qref_cryst.dat in the parent directory.'''

        subdir = self.subdir

        cwd = os.getcwd()
        print "  Calculating native contacts for: ",subdir
        if os.path.exists(cwd+"/"+subdir+"/Qref_cryst.dat"):
            print "  Native contact map "+subdir+"/Qref_cryst.dat exists."
            print "  Skipping shadow map calculation. Loading native contact map..."
            Qref = np.loadtxt(cwd+"/"+subdir+"/Qref_cryst.dat")
            self.contacts = np.loadtxt(cwd+"/"+subdir+"/contacts.dat",dtype=int)
        else:
            print "  Native contact map "+subdir+"/Qref_cryst.dat does not exist."
            print "  Doing shadow map calculation..."
            print "  *** NOTE: module load jdk/1.7.0.21 required for shadow map ***"
            #if not os.path.exists("Qref_shadow"):
            os.chdir(cwd+"/"+subdir+"/Qref_shadow")
            cmd0 = 'cp /projects/cecilia/SCM.1.31.jar .'
            sb.call(cmd0,shell=True,stdout=open("contacts.out","w"),stderr=open("contacts.err","w"))
            #cmd1 = 'echo -e "9\\n3\\n" | pdb2gmx -f clean.pdb -o %s.gro -p %s.top -ignh' % (subdir,subdir)
            cmd1 = 'echo -e "6\\n6\\n" | pdb2gmx -f clean.pdb -o %s.gro -p %s.top -ignh' % (subdir,subdir)
            sb.call(cmd1,shell=True,stdout=open("convert.out","w"),stderr=open("convert.err","w"))
            cmd2 = 'java -jar SCM.1.31.jar -g %s.gro -t %s.top -o %s.contacts -m shadow --coarse CA' % (subdir,subdir,subdir)
            sb.call(cmd2,shell=True,stdout=open("contacts.out","w"),stderr=open("contacts.err","w"))
            self.contacts = np.loadtxt(subdir+".contacts",dtype=int,usecols=(1,3))

            Qref = np.zeros((self.n_residues,self.n_residues))
            for pair in self.contacts:
                Qref[pair[0]-1,pair[1]-1] = 1

            print "  Native contact map calculated with shadow map. Saving Qref_cryst.dat..."
            np.savetxt("Qref_cryst.dat",Qref,delimiter=" ",fmt="%1d")
            np.savetxt(cwd+"/"+subdir+"/Qref_cryst.dat",Qref,delimiter=" ",fmt="%1d")
            np.savetxt(cwd+"/"+subdir+"/contacts.dat",self.contacts,delimiter=" ",fmt="%4d")

        os.chdir(cwd)
        print "  Length = %d  Number of contacts = %d  Nc/L=%.4f" % (len(Qref),sum(sum(Qref)),float(sum(sum(Qref)))/float(len(Qref)))
        self.Qref = Qref
        self.n_contacts = len(self.contacts)
        self.contacts_ndx = "[ contacts ]\n"
        for i in range(self.n_contacts):
            self.contacts_ndx += "%4d %4d\n" % (self.contacts[i][0],self.contacts[i][1])

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='Build a SMOG model.')
    parser.add_argument('--pdb', type=str, required=True, help='pdb')
    args = parser.parse_args() 
    pdb = args.pdb

    model = SmogCalpha(pdb)
