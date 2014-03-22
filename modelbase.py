import numpy as np
import os
import subprocess as sb

''' 
Model Class

Purpose:
    The Model class contains all information of the particular 
structure-based model being used. At this level the Model class is still
general. Each particular model (e.g. Homogeneous Go model, Heterogeneous
Go model, DMC model) has a specific procedure for creating the Gromacs
files so each particular model is a subclass of Model 
(e.g. HomogeneousGoModel).
    This is designed this way so that different models can be easily 
interchanged. Since the model specifics are primarily used to generate
the Gromacs input files (e.g. nonbond_params.itp), the other details
of the model are accessed rarely but are included for completeness and
in the hope that the models will become self-documenting.
    New Model subclasses can be easily written by using the existing 
one as a template. A lot of the functions will be redundant but that is
ok because the most important part is that it is readable.

Description:
    Model contains all the details of the theoretical model we want to 
simulate, like the form and constants of the Hamiltonian. The methods of 
the different Model subclasses are not the most elegant, but they're 
written to maximize readable.

To Do:
- Code a documentation method. Possibly a __repr__ or __str__ method that
  returns a complete summary of the model details in a format that can be
  easily read in later as well.
- Probably going to pop the particular model codes into their own classes.
  (E.g. HomogeneousGoModel.py, DMCModel.py)

'''

class CalphaBase(object):

    def __init__(self):
        pass

    def clean_pdb(self,pdb):
        ''' Clean the given pdb of everything but ATOM lines. Returns an
            all-atom pdb and a C-alpha only pdb.
        PDB fixed-width column format is given by:
        ATOM     44  C   ALA A  11      12.266  21.667  20.517  1.00 28.80           C  
        '''
        first_full = 0
        atomid_full = 1
        cleanpdb_full = ''
        first_ca = 0
        atomid_ca = 1
        cleanpdb_ca = ''
        for line in open(pdb,'r'):
            if line[:3] in ['TER','END']:
                break
            else:
                if line[:4] == 'ATOM':
                    if line[13:15] == "CA":
                        if first_ca == 0:
                            if line[16] in ["A"," "]:
                                newline_ca = 'ATOM%7s  %-4s%3s A%4d%s\n' % \
                                        (atomid_ca,line[13:16],line[17:20],1,line[26:55])
                                atomid_ca += 1
                                first_ca = 1
                                first_index_ca = int(line[22:26]) - 1
                                cleanpdb_ca += newline_ca
                        else:
                            if line[16] in ["A"," "]:
                                newline_ca = 'ATOM%7s  %-4s%3s A%4d%s\n' % \
                                        (atomid_ca,line[13:16],line[17:20],int(line[22:26])-first_index_ca,line[26:55])
                                atomid_ca += 1
                                cleanpdb_ca += newline_ca

                    if first_full == 0:
                        if line[16] in ["A"," "]:
                            newline_full = 'ATOM%7s  %-4s%3s A%4d%s\n' % \
                                    (atomid_full,line[13:16],line[17:20],1,line[26:55])
                            atomid_full += 1
                            first_full = 1
                            first_index_full = int(line[22:26]) - 1
                            cleanpdb_full += newline_full
                    else:
                        if line[16] in ["A"," "]:
                            newline_full = 'ATOM%7s  %-4s%3s A%4d%s\n' % \
                                    (atomid_full,line[13:16],line[17:20],int(line[22:26])-first_index_full,line[26:55])
                            atomid_full += 1
                            cleanpdb_full += newline_full
        
        cleanpdb_full += 'END\n'
        cleanpdb_ca += 'END\n'
        self.cleanpdb = cleanpdb_ca
        self.cleanpdb_full = cleanpdb_full

    def dissect_clean_pdb(self,subdir):
        ''' Extract info from the Native.pdb for making index and 
            itp files '''
        indices = []
        atoms = []
        residues = []
        coords = []

        for line in open(subdir+"/Native.pdb","r"):
            if line.startswith("END"):
                break
            else:
                indices.append(int(line[6:13]))
                atoms.append(line[11:16].strip())
                residues.append(line[17:20])
                coords.append([float(line[31:39]),float(line[39:47]),float(line[47:55])]) 

        ## Coordinates in pdb files are Angstroms. Convert to nanometers.
        coords = np.array(coords)/10.
        return indices, atoms, residues, coords

    def new_get_index_string(self,indices):
        ''' Generates'''
        #ca_indices = atom_indices[j]["CA"]
        ca_string = ''
        i = 1
        for indx in indices: 
            if (i % 15) == 0:
                ca_string += '%4d \n' % indx
            else:
                ca_string += '%4d ' % indx
            i += 1
        ca_string += '\n'
        indexstring = '[ System ]\n'
        indexstring += ca_string
        indexstring += '[ Protein ]\n'
        indexstring += ca_string
        indexstring += '[ Protein-H ]\n'
        indexstring += ca_string
        indexstring += '[ C-alpha ]\n'
        indexstring += ca_string
        indexstring += '[ Backbone ]\n'
        indexstring += ca_string
        indexstring += '[ MainChain ]\n'
        indexstring += ca_string
        indexstring += '[ MainChain+Cb ]\n'
        indexstring += ca_string
        indexstring += '[ MainChain+H ]\n'
        indexstring += ca_string
        indexstring += '[ SideChain ]\n\n'
        indexstring += '[ SideChain-H ]\n\n'
        return indexstring

    def new_get_bonds_itp(self,indices,coords):
        ''' Generate the bonds.itp string.'''
        kb = self.backbone_param_vals["Kb"]
        bonds_string = '[ bonds ]\n'
        for j in range(len(indices)-1):
            i_idx = indices[j]-1
            j_idx = indices[j+1]-1
            dist = np.linalg.norm(coords[i_idx] - coords[j_idx])
            bonds_string += "%5d%5d%5d%16.8f%16.8f\n" %  \
                          (i_idx+1,j_idx+1,1,dist,kb)
        return bonds_string

    def new_get_angles_itp(self,indices,coords):
        ''' Generate the angles.itp string.'''
        ka = self.backbone_param_vals["Ka"]
        angles_string = '[ angles ]\n'
        for j in range(len(indices)-2):
            i_idx = indices[j]-1
            j_idx = indices[j+1]-1
            k_idx = indices[j+2]-1
            xkj = coords[k_idx] - coords[j_idx]
            xkj /= np.linalg.norm(xkj)
            xij = coords[i_idx] - coords[j_idx]
            xij /= np.linalg.norm(xij)
            theta = (180./np.pi)*np.arccos(np.dot(xkj, xij))
            angles_string += "%5d%5d%5d%5d%16.8f%16.8f\n" %  \
                          (i_idx+1,j_idx+1,k_idx+1,1,theta,ka)
        return angles_string

    def new_get_dihedrals_itp(self,indices,coords):
        ''' Write the dihedrals.itp string.'''
        kd = self.backbone_param_vals["Kd"]
        dihedrals_string = '[ dihedrals ]\n'
        dihedrals_ndx_string = '[ dihedrals ]\n'
        for j in range(len(indices)-3):
            i_idx = indices[j]-1
            j_idx = indices[j+1]-1
            k_idx = indices[j+2]-1
            l_idx = indices[j+3]-1
            phi = self.dihedral(coords,i_idx,j_idx,k_idx,l_idx)
            dihedrals_string += "%5d%5d%5d%5d%5d%16.8f%16.8f%5d\n" %  \
                          (i_idx+1,j_idx+1,k_idx+1,l_idx+1,1,phi,kd,1)
            dihedrals_string += "%5d%5d%5d%5d%5d%16.8f%16.8f%5d\n" %  \
                          (i_idx+1,j_idx+1,k_idx+1,l_idx+1,1,3.*phi,kd/2.,3)
            dihedrals_ndx_string += '%4d %4d %4d %4d\n' % (i_idx+1,j_idx+1,k_idx+1,l_idx+1,)
        return dihedrals_string, dihedrals_ndx_string

    def new_get_bonded_itp_strings(self,indices,atoms,residues,coords):
        ''' Create a dictionary of simulation files concerning only 
            bonded interactions. '''
        topology_files = {}
        print "    Creating topol.top"
        topol_top = self.get_topology_string()
        print "    Creating protein.itp"
        protein_itp = self.get_protein_itp()
        print "    Creating bonds.itp"
        bonds_itp = self.new_get_bonds_itp(indices,coords)
        print "    Creating angles.itp"
        angles_itp = self.new_get_angles_itp(indices,coords)
        print "    Creating dihedrals.itp, dihedrals.ndx"
        dihedrals_itp,dihedrals_ndx = self.new_get_dihedrals_itp(indices,coords)
        print "    Creating index.ndx"
        index_ndx = self.new_get_index_string(indices)
        topology_files = {"index.ndx":index_ndx,
                         "topol.top":topol_top,
                         "protein.itp":protein_itp,
                         "bonds.itp":bonds_itp,
                         "angles.itp":angles_itp,
                         "dihedrals.itp":dihedrals_itp,
                         "dihedrals.ndx":dihedrals_ndx}
        return topology_files

    def get_index_string(self,atom_indices):
        ''' Generates'''
        indexs = []
        for j in range(len(atom_indices)):
            ca_indices = atom_indices[j]["CA"]
            ca_string = ''
            i = 1
            for indx in ca_indices: 
                if (i % 15) == 0:
                    ca_string += '%4d \n' % indx
                else:
                    ca_string += '%4d ' % indx
                i += 1
            ca_string += '\n'
            indexstring = '[ System ]\n'
            indexstring += ca_string
            indexstring += '[ Protein ]\n'
            indexstring += ca_string
            indexstring += '[ Protein-H ]\n'
            indexstring += ca_string
            indexstring += '[ C-alpha ]\n'
            indexstring += ca_string
            indexstring += '[ Backbone ]\n'
            indexstring += ca_string
            indexstring += '[ MainChain ]\n'
            indexstring += ca_string
            indexstring += '[ MainChain+Cb ]\n'
            indexstring += ca_string
            indexstring += '[ MainChain+H ]\n'
            indexstring += ca_string
            indexstring += '[ SideChain ]\n\n'
            indexstring += '[ SideChain-H ]\n\n'
            indexs.append(indexstring)
        return indexs

    def get_bonds_itp(self,prots_indices,prots_coords):
        ''' Write the bonds.itp string.'''
        kb = self.backbone_param_vals["Kb"]
        bonds_itps = []
        for i in range(len(prots_indices)):
            indices = prots_indices[i]["CA"]
            coords = prots_coords[i]
            bonds_string = '[ bonds ]\n'
            for j in range(len(indices)-1):
                i_idx = indices[j]-1
                j_idx = indices[j+1]-1
                dist = np.linalg.norm(coords[i_idx] - coords[j_idx])
                bonds_string += "%5d%5d%5d%16.8f%16.8f\n" %  \
                              (i_idx+1,j_idx+1,1,dist/10.,kb)
            bonds_itps.append(bonds_string)
        return bonds_itps

    def get_angles_itp(self,prots_indices,prots_coords):
        ''' Write the angles.itp string.'''
        ka = self.backbone_param_vals["Ka"]
        angles_itps = []
        for i in range(len(prots_indices)):
            indices = prots_indices[i]["CA"]
            coords = prots_coords[i]
            angles_string = '[ angles ]\n'
            for j in range(len(indices)-2):
                i_idx = indices[j]-1
                j_idx = indices[j+1]-1
                k_idx = indices[j+2]-1
                xkj = coords[k_idx] - coords[j_idx]
                xkj /= np.linalg.norm(xkj)
                xij = coords[i_idx] - coords[j_idx]
                xij /= np.linalg.norm(xij)
                theta = (180./np.pi)*np.arccos(np.dot(xkj, xij))
                angles_string += "%5d%5d%5d%5d%16.8f%16.8f\n" %  \
                              (i_idx+1,j_idx+1,k_idx+1,1,theta,ka)
            angles_itps.append(angles_string)
        return angles_itps

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

    def get_dihedrals_itp(self,prots_indices,prots_coords):
        ''' Write the dihedrals.itp string.'''
        kd = self.backbone_param_vals["Kd"]
        dihedrals_itps = []
        dihedrals_ndxs = []
        for i in range(len(prots_indices)):
            indices = prots_indices[i]["CA"]
            coords = prots_coords[i]
            dihedrals_string = '[ dihedrals ]\n'
            dihedrals_ndx_string = '[ dihedrals ]\n'
            for j in range(len(indices)-3):
                i_idx = indices[j]-1
                j_idx = indices[j+1]-1
                k_idx = indices[j+2]-1
                l_idx = indices[j+3]-1
                phi = self.dihedral(coords,i_idx,j_idx,k_idx,l_idx)
                dihedrals_string += "%5d%5d%5d%5d%5d%16.8f%16.8f%5d\n" %  \
                              (i_idx+1,j_idx+1,k_idx+1,l_idx+1,1,phi,kd,1)
                dihedrals_string += "%5d%5d%5d%5d%5d%16.8f%16.8f%5d\n" %  \
                              (i_idx+1,j_idx+1,k_idx+1,l_idx+1,1,3.*phi,kd/2.,3)
                dihedrals_ndx_string += '%4d %4d %4d %4d\n' % (i_idx+1,j_idx+1,k_idx+1,l_idx+1,)
            dihedrals_itps.append(dihedrals_string)
            dihedrals_ndxs.append(dihedrals_ndx_string)
        return dihedrals_itps, dihedrals_ndxs


    def get_itp_strings(self,prots_indices,prots_residues,prots_coords,prots_ndxs,prots_Qref,R_CD=None):
        ''' Create a dictionary of all files needed to run the simulation. 
            These files encompass all the information about the combination
            of the model and the system. The files are returned as a list of
            dictionaries, one dictionary for each each protein in the system.
            Each dictionary holds the following files: *.itp, index.ndx, 
            dihedral.ndx, topol.top, BeadBead.dat. '''
        ## Bonded itp files
        topol_top = self.get_topology_string()
        protein_itp = self.get_protein_itp()
        bonds_itps = self.get_bonds_itp(prots_indices,prots_coords)
        angles_itps = self.get_angles_itp(prots_indices,prots_coords)
        dihedrals_itps,dihedrals_ndxs = self.get_dihedrals_itp(prots_indices,prots_coords)

        ## Move nonbonded file generation to subclass
        atomtypes_itps,atoms_itps = self.get_atomtypes_string(prots_indices,prots_residues)
        nonbond_params_itps,beadbead_files = self.get_nonbond_params_itp(prots_indices,prots_residues,prots_coords,prots_Qref,R_CD=R_CD)

        topology_files = []
        for i in range(len(prots_residues)):
            topology_files.append({"index.ndx":prots_ndxs[i],
                             "topol.top":topol_top,
                             "protein.itp":protein_itp,
                             "atomtypes.itp":atomtypes_itps[i],
                             "atoms.itp":atoms_itps[i],
                             "bonds.itp":bonds_itps[i],
                             "angles.itp":angles_itps[i],
                             "dihedrals.itp":dihedrals_itps[i],
                             "dihedrals.ndx":dihedrals_ndxs[i],
                             "BeadBead.dat":beadbead_files[i],
                             "nonbond_params.itp":nonbond_params_itps[i]})
        return topology_files

    def get_protein_itp(self):
        protein_itp_string = '; molecular topology file for coarse-grained protein\n\n'
        protein_itp_string += '[ moleculetype ]\n'
        protein_itp_string += '; name number.of.exclusions\n'
        protein_itp_string += 'HomGo     3\n\n'
        protein_itp_string += '#include "atoms.itp"\n'
        protein_itp_string += '#include "bonds.itp"\n'
        protein_itp_string += '#include "angles.itp"\n'
        protein_itp_string += '#include "dihedrals.itp"\n'
        return protein_itp_string

    def get_topology_string(self):
        ''' String for topol.top. Need blank line between sections.'''
        topstring = '[ defaults ]\n'
        topstring += '; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n'
        topstring += '  1             1               no              1.0     1.0\n\n'
        topstring += '#include "atomtypes.itp"\n\n'
        topstring += '#include "nonbond_params.itp"\n\n'
        topstring += '#include "protein.itp"\n\n'
        topstring += '[ system ]\n'
        topstring += 'coarse-grained protein\n\n'
        topstring += '[ molecules ]\n'
        topstring +=  'HomGo     1\n'
        return topstring

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
                    "2006, 363, 297-308."}
        return citations[key]

    def residue_mass(self):
        '''Masses in atomic mass units.'''
        residue_mass = {'ALA':   89.0935, 'ARG':  174.2017, 'ASN':  132.1184,
                        'ASP':  133.1032, 'CYS':  121.1590, 'GLN':  146.1451,
                        'GLU':  147.1299, 'GLY':   75.0669, 'HIS':  155.1552,
                        'ILE':  131.1736, 'LEU':  131.1736, 'LYS':  146.1882,
                        'MET':  149.2124, 'PHE':  165.1900, 'PRO':  115.1310,
                        'SER':  105.0930, 'THR':  119.1197, 'TRP':  204.2262,
                        'TYR':  181.1894, 'VAL':  117.1469, 'SOL':   18.0150}
        return residue_mass

    def residue_one_letter_code(self):
        '''Converting from three letter code to one letter FASTA code.'''
        residue_code = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N',
                        'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
                        'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
                        'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                        'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                        'SER': 'S', 'THR': 'T', 'TRP': 'W',
                        'TYR': 'Y', 'VAL': 'V'}
        return residue_code

    def shadow_contacts(self,subdir):
        ''' Call SMOG Shadow jar code to determine the shadow contacts. If 
            the reference matrix Qref_cryst.dat doesn't exist then create 
            and dive into a subdirectory to run shadow map. Then save 
            Qref_cryst.dat in the parent directory.'''

        cwd = os.getcwd()
        print "  Calculating native contacts for: ",subdir
        if os.path.exists(cwd+"/"+subdir+"/Qref_cryst.dat"):
            print "  Native contact map "+subdir+"/Qref_cryst.dat exists."
            print "  Skipping shadow map calculation. Loading native contact map..."
            Qref = np.loadtxt(cwd+"/"+subdir+"/Qref_cryst.dat")
        else:
            print "  Native contact map "+subdir+"/Qref_cryst.dat does not exist."
            print "  Doing shadow map calculation..."
            print "  *** NOTE: module load jdk/1.7.0.21 required for shadow map ***"
            os.chdir(cwd+"/"+subdir+"/Qref_shadow")
            #loadjava = 'module load jdk/1.7.0.21'
            #sb.call(loadjava,shell=True)
            #cmd0 = 'cp /projects/cecilia/ajk8/model_builder/SCM.1.31.jar .'
            cmd0 = 'cp /projects/cecilia/SCM.1.31.jar .'
            sb.call(cmd0,shell=True,stdout=open("contacts.out","w"),stderr=open("contacts.err","w"))
            cmd1 = 'echo -e "9\\n3\\n" | pdb2gmx -f clean.pdb -o %s.gro -p %s.top' % (subdir,subdir)
            sb.call(cmd1,shell=True,stdout=open("convert.out","w"),stderr=open("convert.err","w"))
            cmd2 = 'java -jar SCM.1.31.jar -g %s.gro -t %s.top -o %s.contacts -m shadow --coarse CA' % (subdir,subdir,subdir)
            sb.call(cmd2,shell=True,stdout=open("contacts.out","w"),stderr=open("contacts.err","w"))

            conts = np.loadtxt(subdir+".contacts",usecols=(1,3))
            N = max(conts.ravel())
            Qref = np.zeros((N,N))
            for pair in conts:
                Qref[pair[0]-1,pair[1]-1] = 1

            print "  Native contact map calculated with shadow map. Saving Qref_cryst.dat..."
            np.savetxt("Qref_cryst.dat",Qref,delimiter=" ",fmt="%1d")
            np.savetxt(cwd+"/"+subdir+"/Qref_cryst.dat",Qref,delimiter=" ",fmt="%1d")
        os.chdir(cwd)
        print "  Length = %d  Number of contacts = %d  Nc/L=%.4f" % (len(Qref),sum(sum(Qref)),float(sum(sum(Qref)))/float(len(Qref)))
        self.Qref = Qref

    def write_info_file(self,sub):
        ''' Writes model.info file in subdirectory. The data of the Model object    
            is static, so this is straightforward documentation.'''
        open(sub+"/model.info","w").write(self.__repr__())
