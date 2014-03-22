import numpy as np

from HomogeneousGoModel import HomogeneousGoModel

'''
Mon Mar 3 2014
Alexander Kluber

Heterogeneous Go Model

Purpose:
    The heterogeneous (C-alpha) Go model is very similar to the homogeneous
one, so it is just a sub-class of HomogeneousGoModel with a couple of functions
redefined. 
    It will just take.

Description:

'''

class HeterogeneousGoModel(HomogeneousGoModel):
    ''' All that needs to be computed is the interaction matrix between residues
        i and j.'''

    def __init__(self,contact_energies,disulfides=None,nonbond_param=1.,R_CD=None,dryrun=False):
        self.model_parameters(nonbond_param=nonbond_param,R_CD=R_CD,)
        self.get_interaction_tables()
        self.disulfides = disulfides
        self.cont_type = contact_energies
        self.dryrun = dryrun

    def model_parameters(self,nonbond_param=1.,R_CD=None):
        ''' Contains all the parameter information about the model, as well as
            what types of interactions are included.'''
        self.modelname = "Heterogeneous Go Model"
        self.modelnameshort = "HetGo"
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
        self.R_CD = None
        self.citation = self.citation_info(self.modelnameshort)
    
    def MJ_weights(resi,resj):
        pass

    def Bach_weights(resi,resj):
        pass

    def get_contact_strengths(self,option,path=''):
        ''' Load in the interaction strengths for the desired option. Native contact
            strengths can be taken from MJ parameters, bach parameters, or a saved
            beadbead.dat. Assigns contacts into a Beadbead.dat format i.e. enumerate
            all contact strengths

            

        '''

        if option == "MJ":
            pass
        elif option == "bach":
            pass
        elif option == "load":
            if path != "":
                pass
            else:
                pass

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
                    print resi,resj
                    raise SystemExit
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
