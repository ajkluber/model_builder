""" C-alpha Go-model with non-uniform contact energies.

Description:
    COMING SOON!

Classes:

HomogeneousGoModel
    Variant of homogeneous Go-model (see (1)) with non-uniform contact strengths


References:

(1) Clementi, C.; Nymeyer, H.; Onuchic, J. N. Topological and Energetic
Factors: What Determines the Structural Details of the Transition State
Ensemble and "En-Route" Intermediates for Protein Folding? An Investigation for
Small Globular Proteins. J. Mol. Biol. 2000, 298, 937-953
"""

import numpy as np
import subprocess as sb
import shutil
import os

import mdtraj as md

from HomogeneousGoModel import HomogeneousGoModel

class SecondaryTertiaryGoModel(HomogeneousGoModel):
    """ C-alpha Go-model with non-uniform contact strengths.

    Description:
        
        A subclass of HomogeneousGoModel that treats


    References:

    (1) Clementi, C.; Nymeyer, H.; Onuchic, J. N. Topological and Energetic
    Factors: What Determines the Structural Details of the Transition State
    Ensemble and "En-Route" Intermediates for Protein Folding? An Investigation for
    Small Globular Proteins. J. Mol. Biol. 2000, 298, 937-953

    """

    def __init__(
            self,contact_energies,disulfides=None,nonbond_param=1.,
            cutoff=None,R_CD=None,epsilon_bar=None,dryrun=False):
        self.model_parameters(nonbond_param=nonbond_param,R_CD=R_CD,epsilon_bar=epsilon_bar)
        self.get_interaction_tables()
        self.disulfides = disulfides
        self.cutoff = cutoff
        self.contact_energies = contact_energies
        self.dryrun = dryrun

    def model_parameters(self,nonbond_param=1.,R_CD=None,epsilon_bar=None):
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
        self.epsilon_bar = epsilon_bar
        self.citation = self.citation_info(self.modelnameshort)

    def extract_backbone_Hbonds(self,System):
        """ Use Gromacs and MDTraj to get backbone hydrogen bonds """

        name = System.subdir

        cwd = os.getcwd()
        
        if not os.path.exists(System.path+"/"+name+"/hbonds"):
            os.makedirs(System.path+"/"+name+"/hbonds")
        if not os.path.exists(System.path+"/"+name+"/hbonds/hbonds.dat"):
            ## Determine backbone H-bonds with Gromacs and MDTraj
            shutil.copy(name+".pdb",System.path+"/"+name+"/hbonds/")
            os.chdir(System.path+"/"+name+"/hbonds")
            prep_script = "#!/bin/bash\n"
            prep_script += "echo -e '1\\n6\\n' | pdb2gmx -f %s.pdb -o %s.gro\n" % (name,name)
            prep_script += "trjconv -f %s.gro -o %s.xtc\n" % (name,name)
            prep_script += "editconf -f %s.gro -o %s_withH.pdb\n" % (name,name)
            open("prep_pdb.sh","w").write(prep_script)
            sb.call("bash prep_pdb.sh", shell=True)

            traj = md.load("%s.xtc" % name, top="%s_withH.pdb" % name)
            M = (md.kabsch_sander(traj))[0]
            Hbonds = np.array(M.todense())
            np.savetxt("hbonds.dat",Hbonds)
            os.chdir("..")
        else:
            M = np.loadtxt("hbonds/hbonds.dat")
        os.chdir(cwd)

        self.Hbonds = Hbonds

    def get_contact_strengths(self,indices,residues):
        """ Calculate contact strengths depending on secondary + tertiary contacts
        
        Description:
    
            Add additional contacts depending on if residues i and j are H-bonding.
        Increase the strength 

        """

        N = len(self.Qref)
        contact_epsilons = np.zeros(int(((N-2)*(N-3))/2),float)
        k = 0
        for i in range(len(residues)):
            for j in range(i+3,len(residues)):
                eps = 0
                if self.Hbonds[i][j] < -0.5:
                    eps += 0.3
                if self.Hbonds[j][i] < -0.5:
                    eps += 0.3
                if self.Hbonds[j][i-1] < -0.5:
                    eps += 0.3
                if self.Hbonds[j][i+1] < -0.5:
                    eps += 0.3
                if self.Qref[i][j] == 1:
                    eps += 1.0
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

            for j in range(i+3,len(indices)):
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
                elif contact_epsilons[k] != 0:
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
                    c12 = 5.0*(sig**12)
                    c10 = 6.0*(sig**10)*delta
                    interaction_num = '0'
                k += 1
                beadbead_string += '%5d%5d%8s%8s%5s%16.8E%16.8E%16.8E\n' % \
                        (i_idx,j_idx,resi+str(i_idx),resj+str(j_idx),interaction_num,sig,Knb,delta)
                nonbond_params_string += "%8s%8s%3d  %10.8e  %10.8e\n" % \
                        (resi+str(i_idx),resj+str(j_idx),1,c10,c12)
        return nonbond_params_string,beadbead_string

    
    def prepare_system(self,System):
        ''' Extract all the topology files from Model. '''

        print "Preparing input files for subdirectory:", System.subdir
        print "  Extracting backbone H-bonds..."
        self.extract_backbone_Hbonds(System)

        print "  Cleaning pdb..."
        self.clean_pdb(System.path+"/"+System.subdir+".pdb")
        open(System.path+"/"+System.subdir+"/Native.pdb","w").write(self.cleanpdb)
        open(System.path+"/"+System.subdir+"/Qref_shadow/clean.pdb","w").write(self.cleanpdb_full)
        open(System.path+"/"+System.subdir+"/clean.pdb","w").write(self.cleanpdb_full)
        open(System.path+"/"+System.subdir+"/clean_noH.pdb","w").write(self.cleanpdb_full_noH)
        self.shadow_contacts(System.subdir)

        ## Get index files. Should all be done by modelbase/Calphabase class.
        print "  Dissecting Native.pdb for: indices, atom names, residue names, and coordinates..."
        indices,atoms,residues,coords = self.dissect_clean_pdb(System.subdir)
        print "  Generating bonded files:"
        topology_files = self.get_bonded_itp_strings(indices,atoms,residues,coords)

        #print self.n_contacts,self.n_residues  ## DEBUGGING

        ## Make sure specified disulfides are reasonable: close in space and corresponding with 
        ## cysteine pair.
        self.check_disulfides(residues,coords)

        print "  Generating nonbonded files:"
        atomtypes_itp, atoms_itp = self.get_atomtypes_string(indices,residues)
        nonbond_params_itp, BeadBead_dat = self.get_nonbonded_itp_strings(indices,atoms,residues,coords)

        topology_files["Native.pdb"] = self.cleanpdb
        topology_files["atomtypes.itp"] = atomtypes_itp
        topology_files["atoms.itp"] = atoms_itp
        topology_files["BeadBead.dat"] = BeadBead_dat
        topology_files["nonbond_params.itp"] = nonbond_params_itp
        print "Done preparing:", System.subdir

        System.topology_files = topology_files
