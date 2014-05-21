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
    
    def get_nonbonded_itp_strings(self,indices,atoms,residues,coords):
        """ Get the nonbond_params.itp and BeadBead.dat strings. """
    
        """
            Differentiate between secondary structure interactions (i.e.
        residues interacting by backbone hydrogen-bonds) from tertiary
        interactions (i.e. residues interacting by sidechain contacts).

        """
        pass
        
    
