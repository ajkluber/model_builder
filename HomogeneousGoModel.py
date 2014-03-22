import numpy as np

from modelbase import CalphaBase 


'''
Nov 2013
Alexander Kluber

Homogeneous Go Model

Purpose:
    This class generates the non-bonded interactions for the homogeneous go
model.  It inherihits a lot of functionality from CalphaBase. In this model the
only non-bonded interactions are attractive native interactions.

Description:

'''

class HomogeneousGoModel(CalphaBase):
    ''' This model is the classic Go-model used by Clementi, et.al.(2000).
        Attractive LJ12-10 native contacts and repulsive LJ12 non-contacts.
    '''

    def __init__(self,disulfides=None,nonbond_param=1.,R_CD=None,cutoff=None,dryrun=False):
        self.model_parameters(nonbond_param=nonbond_param,R_CD=R_CD)
        self.get_interaction_tables()
        self.disulfides = disulfides
        self.cutoff = cutoff
        self.contact_energies = None
        self.dryrun = dryrun

    def __repr__(self):
        ''' The string representation of all the model info.'''
        repstring = "[ Model_Name ]\n"
        repstring += "%s\n" % self.modelname
        repstring += "[ Model_Code ]\n" 
        repstring += "%s\n" % self.modelnameshort
        repstring += "[ Bead_Model ]\n" 
        repstring += "%s\n" % self.beadmodel[0]
        repstring += "[ Interaction_Groups ]\n" 
        repstring += "%s\n" % self.interaction_groups[0]
        repstring += "[ Interaction_Types ]\n" 
        repstring += "%s\n" % self.interaction_types[0]
        repstring += "[ Backbone_params ]\n" 
        repstring += "  %5s       %5s       %5s\n" % ("Kb","Ka","Kd")
        repstring += "[ Backbone_param_vals ]\n" 
        repstring += "%10.2f%10.2f%10.2f\n" % (self.backbone_param_vals["Kb"],self.backbone_param_vals["Ka"],self.backbone_param_vals["Kd"])
        repstring += "[ nonbond_param ]\n" 
        repstring += "%10f\n" % self.nonbond_param
        repstring += "[ Solvent ]\n"
        repstring += "%s\n" % self.solvent
        repstring += "[ Disulfides ]\n"
        if self.disulfides == None:
            repstring += "%s\n" % None
        else:
            temp = ''
            for x in self.disulfides:
                temp += " %d " % x 
            repstring += "%s\n" % temp
        repstring += "[ R_CD ]\n"
        repstring += "%s\n" % str(self.R_CD)
        repstring += "[ Cutoff ]\n"
        repstring += "%s\n" % str(self.cutoff)
        repstring += "[ Contact_Energies ]\n"
        repstring += "%s\n" % str(self.contact_energies)
        repstring += "[ Reference ]\n" 
        repstring += "%s\n" % self.citation
        return repstring

    def model_parameters(self,nonbond_param=1.,R_CD=None):
        ''' Contains all the parameter information about the model, as well as
            what types of interactions are included.'''
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
        self.nonbond_param = nonbond_param
        self.R_CD = R_CD
        self.citation = self.citation_info(self.modelnameshort)

    def get_interaction_tables(self):
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


    def get_table_Protein_Protein(self,r):
        ''' Protein-Protein interactions with LJ12-10''' 
        table = np.zeros((len(r),6),float)
        table[:,0] = 1./r
        table[:,1] = 1./(r**2)
        table[:,2] = -1./(r**10)
        table[:,3] = -10./(r**11)
        table[:,4] = 1./(r**12)
        table[:,5] = 12./(r**13)
        return table

    def new_get_atomtypes_string(self,indices,residues):
        ''' Generate atomtypes string.'''
        print "    Creating atomtypes.itp, atoms.itp"
        atoms_itp = '[ atoms ]\n'
        atoms_itp += ' ;  nr  type  resnr  residue  atom  cgnr   charge  mass \n'
        atomtypes_itp = '[ atomtypes ]\n'
        atomtypes_itp += ' ;  name  index      mass         charge  ptype     c10           c12\n'
        for j in range(len(indices)):
            resname = residues[j]
            ## Potentially just one one-letter code here for the atomtypes.
            ## rescode = residue_code[resname]
            index = str(indices[j])
            mass = 100.0
            charge = 0.0 
            ptype = "A"
            atomtypes_itp += "%8s%5s    %10.6f    %10.6f  %s    %10.6f    %10.6f\n" % \
                    (resname+index,index,mass,charge,ptype,0.,0.)
            atoms_itp += "%5s%8s%5s%5s%8s%5s    %10.6f    %10.6f\n" % \
                    (index,resname+index,index,resname,resname+index,index,charge,mass)
        return atomtypes_itp, atoms_itp

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

    def get_nonbond_sigma(self,resi,resj,delta,xi,xj):
        ''' Extract the equilibrium non-bonded interaction distance from the 
            histogram files. NOT USED'''

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

    def new_get_nonbonded_itp_strings(self,indices,atoms,residues,coords):
        ''' Get the nonbond_params.itp and BeadBead.dat strings. '''
    
        ## Options to check:
        ## R_CD
        ## disulfides
        
        print "    Creating nonbond_params.itp, BeadBead.dat"

        if self.R_CD != None:
            print "    Using R_C/D option: ", self.R_CD
            Nc = float(sum(sum(self.Qref)))
            Nd = float(len(self.Qref)-4)
            Knb = (self.R_CD*Nd/Nc)*self.backbone_param_vals["Kd"]
            self.nonbond_param = Knb 
        else:
            Knb = self.nonbond_param

        print "    Nonbonded multiplier: ", Knb
        
        print "    Disulfides: ", self.disulfides
        native = 0
        interaction_counter = 1
        nonbond_params_string = '[ nonbond_params ]\n'
        beadbead_string = ''
        for i in range(len(indices)):
            if self.disulfides != None:
                if (i+1) in self.disulfides[::2]:
                    #print self.disulfides[::2]     ## DEBUGGING
                    ds_flag = 1
                    partner = self.disulfides[self.disulfides.index(i+1) + 1] - 1
                    #print "## Residue ", i+1, " is in a disulfide with ", partner+1    ## DEBUGGING
                else:
                    ds_flag = 0
            else:
                ds_flag = 0

            for j in range(i+4,len(indices)):
                resi = residues[i]
                resj = residues[j]
                i_idx = indices[i]
                j_idx = indices[j]
                delta = j - i
                xi = coords[i]
                xj = coords[j]
                if (ds_flag == 1) and (j == partner):
                    ## Making the link between disulfides stronger. Special interaction number.
                    print "    Linking disulfides between residues: ",residues[i]+str(i+1)," and ",residues[partner]+str(partner+1)
                    c12 = self.backbone_param_vals["Kb"]*5.0*(sig**12)
                    c10 = self.backbone_param_vals["Kb"]*6.0*(sig**10)*delta
                    interaction_num = 'ss'
                elif self.Qref[i][j] == 1:
                    ## Regular native contact are attractive.
                    sig = np.linalg.norm(xi - xj)
                    delta = 1
                    c12 = Knb*5.0*(sig**12)
                    c10 = Knb*6.0*(sig**10)*delta
                    interaction_num = str(interaction_counter)
                    interaction_counter += 1
                else:
                    ## Non-native interactions are repulsive at constant distance of 3.5A.
                    #sig, delta = self.get_nonbond_sigma(resi,resj,delta,xi,xj)  ## DEPRECATED
                    sig = 0.35
                    delta = 0
                    c12 = Knb*5.0*(sig**12)
                    c10 = Knb*6.0*(sig**10)*delta
                    interaction_num = '0'
                native += delta
                beadbead_string += '%5d%5d%8s%8s%5s%16.8E%16.8E%16.8E\n' % \
                        (i_idx,j_idx,resi+str(i_idx),resj+str(j_idx),interaction_num,sig,Knb,delta)
                nonbond_params_string += "%8s%8s%3d  %10.8e  %10.8e\n" % \
                        (resi+str(i_idx),resj+str(j_idx),1,c10,c12)
                         #(resi+str(i_idx),resj+str(j_idx),1,sig,c12) ## DEBUGGING
            #print native   ## DEBUGGING
            #print nonbond_params_string ## DEBUGGING
            #raise SystemExit
        return nonbond_params_string,beadbead_string

    def get_nonbond_params_itp(self,prots_indices,prots_residues,prots_coords,prots_Qref,R_CD=None):
        ''' Get the nonbond_params.itp and BeadBead.dat strings. Select bond
            distances for native contacts. This is the core of what 
            distinquishes the different models. Also collects the native 
            contacts. '''
    
        nonbond_itps = []
        beadbead_files = []
        for prot_num in range(len(prots_indices)):
            indices = prots_indices[prot_num]["CA"]
            residues = prots_residues[prot_num]
            coords = prots_coords[prot_num]/10.

            if R_CD != None:
                Nc = float(sum(sum(prots_Qref[prot_num])))
                Nd = float(len(prots_Qref[prot_num])-4)
                Knb = (R_CD*Nd/Nc)*self.backbone_param_vals["Kd"]
            else:
                Knb = self.nonbond_param
            
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
                    if prots_Qref[prot_num][i][j] == 1:
                        ## Native contact are attractive.
                        sig = np.linalg.norm(xi - xj)
                        delta = 1
                        c12 = Knb*5.0*(sig**12)
                        c10 = Knb*6.0*(sig**10)*delta
                    else:
                        ## Non-native interactions are repulsive at constant distance of 3.5A.
                        #sig, delta = self.get_nonbond_sigma(resi,resj,delta,xi,xj)
                        sig = 0.35
                        delta = 0
                        c12 = Knb*5.0*(sig**12)
                        c10 = Knb*6.0*(sig**10)*delta
                    native += delta
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


    def check_disulfides(self,residues,coords):
        ''' Check that specified disulfides are between cysteine and that 
            the corresonding pairs are within 0.8 nm.'''
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

        else:
            print "  No disulfides to check."

    def new_prepare_system(self,System):
        ''' Extract all the topology files from Model. 
        '''

        print "Preparing input files for subdirectory:", System.subdir
        print "  Cleaning pdb..."
        self.clean_pdb(System.path+"/"+System.subdir+".pdb")
        open(System.path+"/"+System.subdir+"/Native.pdb","w").write(self.cleanpdb)
        open(System.path+"/"+System.subdir+"/Qref_shadow/clean.pdb","w").write(self.cleanpdb_full)
        open(System.path+"/"+System.subdir+"/clean.pdb","w").write(self.cleanpdb_full)
        self.shadow_contacts(System.subdir)

        ## Get index files. Should all be done by modelbase/Calphabase class.
        print "  Dissecting Native.pdb for: indices, atom names, residue names, and coordinates..."
        indices,atoms,residues,coords = self.dissect_clean_pdb(System.subdir)
        print "  Generating bonded files:"
        topology_files = self.new_get_bonded_itp_strings(indices,atoms,residues,coords)

        ## Make sure disulfide
        self.check_disulfides(residues,coords)

        print "  Generating nonbonded files:"
        atomtypes_itp, atoms_itp = self.new_get_atomtypes_string(indices,residues)
        nonbond_params_itp, BeadBead_dat = self.new_get_nonbonded_itp_strings(indices,atoms,residues,coords)

        topology_files["Native.pdb"] = self.cleanpdb
        topology_files["atomtypes.itp"] = atomtypes_itp
        topology_files["atoms.itp"] = atoms_itp
        topology_files["BeadBead.dat"] = BeadBead_dat
        topology_files["nonbond_params.itp"] = nonbond_params_itp

        #print coords  ## DEBUGGING
        #print indices, atoms, residues, coords  ## DEBUGGING
        #print atoms_itp        ## DEBUGGING
        #print atomtypes_itp        ## DEBUGGING
        #print nonbond_params_itp        ## DEBUGGING
        print "Done preparing:", System.subdir

        System.topology_files = topology_files
