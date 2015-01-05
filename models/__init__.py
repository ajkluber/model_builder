""" Submodule with classes of coarse-grain models.

Description:

    A module that prepares the Model object which contains the functions and
parameters needed to prepare input files for coarse-grain protein simulations
in Gromacs. Each model is a seperate class. 


Classes:

    SmogCalpha


References:


"""

import SmogCalpha
import bonded_potentials
import pairwise_potentials
import pdb_parser
import beadbead_to_params
#import prepare_compound_pairwise

def citation_info(key):
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
