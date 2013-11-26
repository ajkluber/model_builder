import numpy as np

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


        
## Some dictionaries of helpful information.
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

## Masses in atomic mass units.
residue_mass = {'ALA':   89.0935, 'ARG':  174.2017, 'ASN':  132.1184,
                'ASP':  133.1032, 'CYS':  121.1590, 'GLN':  146.1451,
                'GLU':  147.1299, 'GLY':   75.0669, 'HIS':  155.1552,
                'ILE':  131.1736, 'LEU':  131.1736, 'LYS':  146.1882,
                'MET':  149.2124, 'PHE':  165.1900, 'PRO':  115.1310,
                'SER':  105.0930, 'THR':  119.1197, 'TRP':  204.2262,
                'TYR':  181.1894, 'VAL':  117.1469, 'SOL':   18.0150}

## Converting from three letter code to one letter FASTA code.
residue_code = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N',
                'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
                'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
                'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W',
                'TYR': 'Y', 'VAL': 'V'}

class HomogeneousGoModel(object):
    ''' This model is the closest model to the classical structure-based model 
        used by Clementi, et.al.(2000). The native contacts are the only 
        attractive terms in the non-bonded potential. The non-native contacts
        are treated with the repulsive r**12 term. This model is not exactly the
        same as the original SBMs in that the non-native interactions still have
        a "shape" in the same way as the other Clementi models.'''

    def __init__(self,path):
        self.path = path
        self.model_parameters()
        self.get_tables()

    def get_tables(self):
        ''' Returns the table for interaction type 'i'. The values of the
            potential is set to 0.0 when r is too close to zero to avoid
            blowup. Writes interaction function table files: 
            table_GROUP1_GROUP2.xvg for each pair [GROUP1,GROUP2] listed in
            self.interaction_groups.'''
        r = np.arange(0.0,4.0,0.002)
        self.tables = []
        for i in range(len(self.interaction_groups)):
            table = np.zeros((len(r),7),float)
            table[:,0] = r
            table[20:,1:7] = self.interaction_tables[i](r[20:])
            self.tables.append(table)
        self.other_table = np.zeros((len(r),7),float)
        self.other_table[:,0] = r

    def __repr__(self):
        ''' The string representation of all the model info.'''
        repstring = "Model Name: %s\n" % self.modelname
        repstring += "Reference: %s\n" % self.citation
        repstring += "Model Code: %s\n" % self.modelnameshort
        repstring += "Bead Model : %s\n" % self.beadmodel[0]
        repstring += "Interaction Groups: %s\n" % self.interaction_groups[0]
        repstring += "Interaction Types: %s\n" % self.interaction_types[0]
        repstring += "Solvent: %s\n" % self.solvent
        return repstring

    def model_parameters(self):
        self.modelname = "Homogenous Go Model"
        self.modelnameshort = "HomGo"
        self.beadmodel = ["CA"]
        self.energygroups = ["Protein"]
        self.energygrps = ""
        for grp in self.energygroups: self.energygrps += grp + " "
        #self.energygrps = ["Protein"]
        self.interaction_groups = ["Protein_Protein"]
        self.interaction_types = ["LJ12-10"]
        self.interaction_tables = [ self.get_table_Protein_Protein ]
        self.energygrp_table = ""
        for intgrps in self.interaction_groups: 
            self.energygrp_table += intgrps.split("_")[0]+" "+intgrps.split("_")[1]+" "
        self.solvent = "None"
        self.backbone_params = ["Kb","Ka","Kd"]
        self.backbone_param_vals = {"Kb":20000.,"Ka":400.,"Kd":1}
        self.nonbond_param = 1.
        self.citation = citations[self.modelnameshort]

    def get_table_Protein_Protein(self,r):
        ''' protein-protein interactions with LJ12-10''' 
        table = np.zeros((len(r),6),float)
        table[:,0] = 1./r
        table[:,1] = 1./(r**2)
        table[:,2] = -1./(r**10)
        table[:,3] = -10./(r**11)
        table[:,4] = 1./(r**12)
        table[:,5] = 12./(r**13)
        return table
        
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

    def get_index_string(self,subdirs,atom_indices):
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

    def get_atomtypes_string(self,atom_indices,residue_names):
        ''' Create atomtypes string.'''
        atomtypes_itps = []
        atoms_itps = []
        for i in range(len(residue_names)):
            atoms = '[ atoms ]\n'
            atoms += ' ;  nr  type  resnr  residue  atom  cgnr   charge  mass \n'
            atomtypes = '[ atomtypes ]\n'
            atomtypes += ' ;  name  index      mass         charge  ptype     c10           c12\n'
            res_names = residue_names[i]
            indices = atom_indices[i]["CA"]
            for j in range(len(res_names)):
                resname = res_names[j]
                ## Potentially just one one-letter code here for the atomtypes.
                ## rescode = residue_code[resname]
                ndx = str(indices[j])
                mass = 100.0
                charge = 0.0 
                ptype = "A"
                atomtypes += "%8s%5s    %10.6f    %10.6f  %s    %10.6f    %10.6f\n" % \
                        (resname+ndx,ndx,mass,charge,ptype,0.,0.)
                atoms += "%5s%8s%5s%5s%8s%5s    %10.6f    %10.6f\n" % \
                        (ndx,resname+ndx,ndx,resname,resname+ndx,ndx,charge,mass)
            atomtypes_itps.append(atomtypes)
            atoms_itps.append(atoms)
        return atomtypes_itps, atoms_itps

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

    def get_nonbond_sigma(self,resi,resj,delta,xi,xj):

        histpath = "/projects/cecilia/ajk8/model_builder/Histograms/contact_CA/hist"
        names = sorted([resi,resj])
        if delta == 4:
            extpath = "4"
        elif delta == 5 or delta == 6:
            extpath = "56" 
        elif delta >= 7:
            extpath = "7up"
        else:
            print "You should never see this"
        histpath += extpath+"/"+names[0]+"_"+names[1]+".xvg"

        x, Px = np.loadtxt(histpath,unpack=True)
        cumx = np.array([ sum(Px[:i]) for i in range(len(Px)) ])
        cumx[-1] = 1
        q = np.random.random()
        index = list(cumx > q).index(True)
        sig = 0.05*(x[index-1] + x[index])
        natsig = np.linalg.norm(xi-xj)
        newsig = min([natsig,sig])
        #print natsig, newsig ## DEBUGGING
        if natsig <= 1.25*newsig:
            native = 1
        else:
            native = 0
        return newsig, native


    def get_nonbond_params_itp(self,prots_indices,prots_residues,prots_coords):
        ''' Get the nonbond_params.itp and BeadBead.dat strings. Select bond
            distances for native contacts. This is the core of what 
            distinquishes the different models. Also collects the native 
            contacts.'''
    
        Knb = self.nonbond_param
        nonbond_itps = []
        beadbead_files = []
        for prot_num in range(len(prots_indices)):
            indices = prots_indices[prot_num]["CA"]
            residues = prots_residues[prot_num]
            coords = prots_coords[prot_num]/10.
    
            beadbead_string = ''
            native = 0
            interaction_num = 1
            nonbond_params_string = '[ nonbond_params ]\n'
            for i in range(len(indices)):
                for j in range(i+4,len(indices)):
                    resi = residues[i]
                    resj = residues[j]
                    i_idx = indices[i]
                    j_idx = indices[j]
                    delta = j - i
                    xi = coords[i]
                    xj = coords[j]
                    sig, delta = self.get_nonbond_sigma(resi,resj,delta,xi,xj)
                    native += delta
                    c12 = Knb*5.0*(sig**12)
                    c10 = Knb*6.0*(sig**10)*delta
                    beadbead_string += '%5d%5d%8s%8s%5d%16.8E%16.8E%16.8E\n' % \
                            (i_idx,j_idx,resi+str(i_idx),resj+str(j_idx),interaction_num,sig,Knb,delta)
                    nonbond_params_string += "%8s%8s%3d  %10.8e  %10.8e\n" % \
                            (resi+str(i_idx),resj+str(j_idx),1,c10,c12)
                             #(resi+str(i_idx),resj+str(j_idx),1,sig,c12) ## DEBUGGING
                    interaction_num += 1
            #print native   ## DEBUGGING
            #print nonbond_params_string ## DEBUGGING
            #raise SystemExit
            nonbond_itps.append(nonbond_params_string)
            beadbead_files.append(beadbead_string)

        return nonbond_itps,beadbead_files

    def get_itp_strings(self,prots_indices,prots_residues,prots_coords,prots_ndxs):
        ''' Create a dictionary of all files needed to run the simulation. 
            These files encompass all the information about the combination
            of the model and the system. The files are returned as a list of
            dictionaries, one dictionary for each each protein in the system.
            Each dictionary holds the following files: *.itp, index.ndx, 
            dihedral.ndx, topol.top, BeadBead.dat. '''
        topol_top = self.get_topology_string()
        protein_itp = self.get_protein_itp()
        atomtypes_itps,atoms_itps = self.get_atomtypes_string(prots_indices,prots_residues)
        bonds_itps = self.get_bonds_itp(prots_indices,prots_coords)
        angles_itps = self.get_angles_itp(prots_indices,prots_coords)
        dihedrals_itps,dihedrals_ndxs = self.get_dihedrals_itp(prots_indices,prots_coords)
        nonbond_params_itps,beadbead_files = self.get_nonbond_params_itp(prots_indices,prots_residues,prots_coords)

        #print bonds_itps[0]     ##DEBUGGING
        #print angles_itps[0]    ##DEBUGGING
        #print dihedrals_itps[0] ##DEBUGGING
        #raise SystemExit
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

    def write_info_file(self):
        ''' Writes model.info file for new simulation project. 
            DOESN"T WORK. NEEDS TO BE UPDATED.'''

        template_info = open('/home/ajk8/projects/dmc_model/gmx/%s.params' % self.modeltype, 'r').read()
        pdbs = ''
        for pdb in self.pdbs: pdbs += pdb + " "
        info = open('model.info','w')
        info.write("Simulation project initialized on: %s\n" % time.asctime())
        info.write("[ System ]\n")
        info.write("SystemName:    %s\n"  % self.systemname )
        info.write("PDBs:    %s\n"  % pdbs )
        info.write("\n")
        info.write("[ Model ]\n")
        info.write(template_info)
        info.write("Solvent: %s" % self.solvent)
        info.close()
