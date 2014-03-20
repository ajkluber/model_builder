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
    ''' This model is the closest model to the classical structure-based model 
        used by Clementi, et.al.(2000). The native contacts are the only 
        attractive terms in the non-bonded potential. The non-native contacts
        are treated with the repulsive r**12 term. This model is not exactly the
        same as the original SBMs in that the non-native interactions still have
        a "shape" in the same way as the other Clementi models.'''

    def __init__(self,disulfides=None,nonbond_param=1.,R_CD=None,cutoff=None):
        self.model_parameters(nonbond_param=nonbond_param,R_CD=R_CD)
        self.get_tables()
        self.disulfides = disulfides
        self.cutoff = cutoff
        self.contact_energies = None

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

    def load_info_file(self,sub):
        ''' NOT IN USE.'''
        
        info_file = open(sub+'/model.info','r')
        line = info_file.readline()
        while line != '':
            #print line[:-1]
            print line.split()
            value = line.split()[1]
            if value == "":
                self.Tf_iteration[i] = int(info_file.readline())
            elif value == "System_Name":
                self.systemname = info_file.readline()[:-1]
            else:
                if line[0] == "[":
                    line = info_file.readline()
            line = info_file.readline()

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


    def prepare_system_new(System,options):
        ''' Extract all the topology files from Model. 
            options:
                disulfides
                R_CD
                cutoff
                dryrun
        '''
        print "Preparing files..."
        prots_Qref = System.shadow_contacts()
        System.write_Native_pdb_CA()
        if R_CD != None:
            for i in range(len(System.subdirs)):
                Nc = float(sum(sum(prots_Qref[i])))
                Nd = float(len(prots_Qref[i])-4)
                System.nonbond_params.append((R_CD*Nd/Nc)*self.backbone_param_vals["Kd"])
                System.R_CD.append(R_CD)
        else:
            for i in range(len(System.subdirs)):
                N = len(prots_Qref[i])
                Nc = float(sum(sum(prots_Qref[i])))
                Nd = float(len(prots_Qref[i])-4)
                print "Num contacts per residue: ",Nc/N
                System.nonbond_params.append(self.nonbond_param)
                System.R_CD.append(None)
        prots_indices, prots_residues, prots_coords = System.get_atom_indices(self.beadmodel)
        prots_ndxs = self.get_index_string(prots_indices)
        topology_files = self.get_itp_strings(prots_indices, prots_residues, prots_coords,prots_ndxs,prots_Qref,R_CD=R_CD)
        System.topology_files = topology_files
