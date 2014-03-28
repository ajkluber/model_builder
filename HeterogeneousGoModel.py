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

    def __init__(self,contact_energies,disulfides=None,nonbond_param=1.,cutoff=None,R_CD=None,dryrun=False):
        self.model_parameters(nonbond_param=nonbond_param,R_CD=R_CD,)
        self.get_interaction_tables()
        self.disulfides = disulfides
        self.cutoff = cutoff
        self.contact_energies = contact_energies
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
    
    def get_MJ_weights(self,resi,resj):
        return 1.

    def get_Bach_weights(self,resi,resj):
        return 1.

    def get_MC2004_weights(self):
        ''' Load in contact strengths from BeadBead.dat file.'''
        beadbead = np.loadtxt(self.contact_energies,dtype=str)
        contact_epsilons = beadbead[:,6].astype(float)
        return contact_epsilons

    def get_contact_strengths(self,indices,residues):
        ''' Load in the interaction strengths for the desired option. Native contact
            strengths can be taken from MJ parameters, bach parameters, or a saved
            beadbead.dat. Assigns contacts into a Beadbead.dat format i.e. enumerate
            all contact strengths

        '''
        if self.contact_energies.endswith("BeadBead.dat"):
            ## Load contact strengths from file.
            contact_epsilons = self.get_MC2004_weights()
        else:
            ## Load contact strengths.
            N = len(self.Qref)
            contact_epsilons = np.zeros(int(((N-3)*(N-4))/2),float)
            k = 0
            for i in range(len(residues)):
                for j in range(i+4,len(residues)):
                    if self.Qref[i][j] == 1:
                        resi = residues[i]
                        resj = residues[j]
                        if self.contact_energies == "MJ":
                            eps = self.get_MJ_weights(resi,resj) 
                        elif self.contact_energies == "Bach":
                            eps = self.get_Bach_weight(resi,resj)
                        elif self.contact_energies == "MC2004":
                            eps = 1.0
                        else:
                            print "ERROR!"
                            print "  Couldn't find contact_energies for:", self.contact_energies
                            print "  Exiting."
                            raise SystemExit
                        contact_epsilons[k] = eps
                    k += 1
        return contact_epsilons
        

    def get_nonbonded_itp_strings(self,indices,atoms,residues,coords):
        ''' Get the nonbond_params.itp and BeadBead.dat strings. '''
    
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

        contact_epsilons = self.get_contact_strengths(indices,residues)
        interaction_counter = 1
        nonbond_params_string = '[ nonbond_params ]\n'
        beadbead_string = ''
        k = 0
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
                    sig = np.linalg.norm(xi - xj)
                    c12 = self.backbone_param_vals["Kb"]*5.0*(sig**12)
                    c10 = self.backbone_param_vals["Kb"]*6.0*(sig**10)*delta
                    interaction_num = 'ss'
                elif self.Qref[i][j] == 1:
                    ## Regular native contact are attractive.
                    
                    Knb = contact_epsilons[k]
                    sig = np.linalg.norm(xi - xj)
                    delta = 1
                    c12 = Knb*5.0*(sig**12)
                    c10 = Knb*6.0*(sig**10)*delta
                    interaction_num = str(interaction_counter)
                    interaction_counter += 1
                else:
                    ## Non-native interactions are repulsive at constant distance of 3.5A.
                    sig = 0.35
                    delta = 0
                    c12 = Knb*5.0*(sig**12)
                    c10 = Knb*6.0*(sig**10)*delta
                    interaction_num = '0'
                k += 1
                beadbead_string += '%5d%5d%8s%8s%5s%16.8E%16.8E%16.8E\n' % \
                        (i_idx,j_idx,resi+str(i_idx),resj+str(j_idx),interaction_num,sig,Knb,delta)
                nonbond_params_string += "%8s%8s%3d  %10.8e  %10.8e\n" % \
                        (resi+str(i_idx),resj+str(j_idx),1,c10,c12)
        return nonbond_params_string,beadbead_string
    
