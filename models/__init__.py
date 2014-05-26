""" Submodule with classes of coarse-grain models.

Description:

    A module that prepares the Model object which contains the functions and
parameters needed to prepare input files for coarse-grain protein simulations
in Gromacs. Each model is a seperate class. 


Classes:



HomogeneousGoModel
    Model with only native contacts and homogeneous contact energies

HeterogeneousGoModel 
    Model with only native contacts and heterogeneous contact energies

DMCModel (coming soon)

References:



"""

import os

import DMCModel
import HeterogeneousGoModel
import HomogeneousGoModel
import SecondaryTertiaryGoModel
import CalphaBase

def check_options(inputoptions):
    ''' Check that all options are compatible and in proper format. Any options
        that are omitted are given the value None (or 1 for nonbond_param) for
        completeness. Returns: options, a dictionary of completed options.'''

    ## List of supported models & corresponding representations.
    available_models = ["HomGo","HetGo","DMC"]
    beadmodels = {"HomGo":["CA"],"HetGo":["CA"]}
    contactopts = {"HetGo":["MJ","Bach","MC2004","FRETFit","SecTer","SecTerSec","SecTerTer"]}

    modelcode = inputoptions["Model_Code"]
    beadmodel = inputoptions["Bead_Model"]
    
    print "Checking that model options are consistent..."
    print "Inputted model options:"
    for key in inputoptions.keys():
        print "  ", key , " = ", inputoptions[key]

    options = {"Model_Code":modelcode, "Bead_Model":beadmodel}
    ## Check if model code is legal.
    if modelcode not in available_models:
        print "ERROR! Invalid Model_Code"
        print "Model: ", modelcode, "  not within available models: ",available_models
        print "Exiting."
        raise SystemExit

    ## Check if bead model is available for model code.
    if beadmodel not in beadmodels[modelcode]:
        print "ERROR! Invalid Bead_Model"
        print "Bead model: ", beadmodel, " not available for model: ", modelcode
        print "Model: ",modelcode," has the following bead models: ", beadmodels[modelcode]
        print "Exiting."
        raise SystemExit

    ## Check for disulfides. If yes disulfids should be in format (e.g.):
    ## [ 1, 6, 7, 20] where the first pair of numbers 1,6 correspond to 
    ## a disulfide bond between those residue indices according to clean.pdb 
    ## indices.
    if inputoptions.has_key("Disulfides"):
        if inputoptions["Disulfides"] not in ["","None",None,False]:
            disulf = inputoptions["Disulfides"]
            if type(disulf) == str:
                disulf = disulf.split()
            if (len(disulf) % 2) != 0:
                print "ERROR! Invalid disulfides argument"
                print "Disulide argument: ",disulf, " must be a list of pairs. Length must be even."
                print "Exiting."
                raise SystemExit
            else:
                try:
                    disulfides = []
                    for dis in disulf:
                        disulfides.append(int(dis))
                except:
                    print "ERROR! Invalid value for disulides"
                    print "Disulfide argument: ",disulf, " must be list of integers"
                    print "Exiting."
                    print SystemExit
        else:
            disulfides = None
    else:
        disulfides = None
    options["Disulfides"] = disulfides

    ## Check for contact energies input. Only valid for certain models. Format
    ## for contact_energies should be an allowed option for the model (e.g. 'MJ'
    ## for 'HetGo') or a path to a Beadbead.dat with previously saved contact
    ## energies.
    if inputoptions.has_key("Contact_Energies"):
        if inputoptions["Contact_Energies"] not in ["","None",None,False]:
            if contactopts.has_key(modelcode):
                if (inputoptions["Contact_Energies"] in contactopts[modelcode]):
                    contact_energies = inputoptions["Contact_Energies"]
                elif inputoptions["Contact_Energies"].endswith("BeadBead.dat"):
                    if not os.path.exists(inputoptions["Contact_Energies"]):
                        print "ERROR!"
                        print "Contact energies option: ", inputoptions["Contact_Energies"], \
                            " points to a nonexistant file!"
                        print "Exiting."
                        raise SystemExit
                    else:
                        contact_energies = inputoptions["Contact_Energies"]
                else:
                    print "ERROR!"
                    print "Model: ",modelcode," contact energies option must be from: ", \
                            contactopts[modelcode], " or in format /path/to/BeadBead.dat"
                    print "Exiting."
                    raise SystemExit
            else:
                print "ERROR!"
                print "Model: ",modelcode," doesn't support contact energies option."
                print "Exiting."
                raise SystemExit
        else:
            if contactopts.has_key(modelcode):
                print "ERROR!"
                print "Model: ",modelcode," requires contact energies option."
                print "Exiting."
                raise SystemExit
            else:
                contact_energies = None
    else:
        contact_energies = None
    options["Contact_Energies"] = contact_energies

    ## Check if average contact strength parameter, epsilon_bar, is set. This 
    ## keeps the average contact strength normalized to some number.
    if inputoptions.has_key("Epsilon_Bar"):
        if inputoptions["Epsilon_Bar"] not in ["","None",None,False]:
            try:
                epsilon_bar = float(inputoptions["Epsilon_Bar"])
            except:
                print "TypeError! epsilon_bar value must be a float!"
                print "Exiting."
                raise SystemExit
        else:
            epsilon_bar = None
    else:
        epsilon_bar = None
    options["Epsilon_Bar"] = epsilon_bar


    ## Check for R_CD option. This option fixes the ratio of contact (C) to
    ## dihedral (D) energy. Soon to be deprecated, replaced by epsilon_bar.
    if inputoptions.has_key("R_CD"):
        if inputoptions["R_CD"] not in ["","None",None,False]:
            try:
                R_CD = float(inputoptions["R_CD"])
            except:
                print "TypeError! R_CD value must be a float!"
                print "Exiting."
                raise SystemExit
        else:
            R_CD = None
    else:
        R_CD = None
    options["R_CD"] = R_CD

    ## This multiplier for the nonbond_param is 1 by default, but may be
    ## different if R_CD was used. This is calculated automatically and not a command
    ## line option.
    if inputoptions.has_key("nonbond_param"):
        if inputoptions["nonbond_param"] not in ["","None",None,False]:
            try:
                nonbond_param = float(inputoptions["nonbond_param"])
            except:
                print "TypeError! nonbond_param value must be a float!"
                print "Exiting."
                raise SystemExit
        else:
            nonbond_param = 1.
    else:
        nonbond_param = 1.
    options["nonbond_param"] = nonbond_param

    ## Cutoff will switch method of determining native contacts from Shadow map to 
    ## using heavy atom contacts within a cutoff of a given residue. Not currently
    ## implemented.
    if inputoptions.has_key("Cutoff"):
        if inputoptions["Cutoff"] not in ["","None",None,False]:
            try:
                cutoff = float(inputoptions["Cutoff"])
            except:
                print "TypeError! cutoff value must be a float!"
                print "Exiting."
                raise SystemExit
        else: 
            cutoff = None
    else:
        cutoff = None
    options["Cutoff"] = cutoff

    ## This option is preemptive in case we want to solvent. Not implemented.
    if inputoptions.has_key("Solvent"):
        if inputoptions["Solvent"] in ["","None",None,False]:
            solvent = None
        else: 
            print "Error! Solvent option not implemented!"
            print "Solvent variable: ", inputoptions["Solvent"], " not allowed. Only: ",["",None]
            print "Exiting."
            raise SystemExit
    else:
        solvent = None
    options["Solvent"] = solvent

    ## Dry run flag will prevent any simulations from being submitted. Used to see if file 
    ## preparation runs smoothly.
    if inputoptions["Dry_Run"] == True:
        dryflag = True
    else:
        dryflag = False
    options["Dry_Run"] = dryflag

    print "Model options cleared!"
    print "Using model options:"
    for key in options.keys():
        print "  ", key , " = ", options[key]
            
    return options

def get_model(options):
    ''' Return a model with the inputted dictionary of options'''
    type = options["Model_Code"]
    if type == "HomGo":
        model = HomogeneousGoModel.HomogeneousGoModel(
                            disulfides=options["Disulfides"],
                            nonbond_param=options["nonbond_param"],
                            R_CD=options["R_CD"],
                            epsilon_bar=options["Epsilon_Bar"],
                            cutoff=options["Cutoff"],
                            dryrun=options["Dry_Run"])
    elif type == "HetGo":
        if options["Contact_Energies"].startswith("SecTer"):
            model = SecondaryTertiaryGoModel.SecondaryTertiaryGoModel(
                            options["Contact_Energies"],
                            disulfides=options["Disulfides"],
                            nonbond_param=options["nonbond_param"],
                            R_CD=options["R_CD"],
                            cutoff=options["Cutoff"],
                            dryrun=options["Dry_Run"])
        else:
            model = HeterogeneousGoModel.HeterogeneousGoModel(
                            options["Contact_Energies"],
                            disulfides=options["Disulfides"],
                            nonbond_param=options["nonbond_param"],
                            R_CD=options["R_CD"],
                            cutoff=options["Cutoff"],
                            dryrun=options["Dry_Run"])
    elif type == "DMC":
        model = DMCModel.DMCModel()
    else:
        ## Due to previous error checking this should never be seen.
        print "ERROR! "
        print "Model: ",type," doesn't exist!"
    return model


def load_model(subdir,dryrun=False):
    ''' Read model.info files in subdirectories and create models.'''
    info_file = open(subdir+'/model.info','r')
    line = info_file.readline()
    options = {"Dry_Run":dryrun}
    while line != '':
        field = line.split()[1]
        value = info_file.readline()
        if field == "Reference":
            break
        elif field in ["Interaction_Groups","Model_Name",
                        "Backbone_params","Backbone_param_vals",
                        "Interaction_Types"]:
            pass
        else:
            options[field] = value[:-1]
        line = info_file.readline()
    options = check_options(options)
    Model = get_model(options)
    return Model

def load_models(subdirs,dryrun=False):
    ''' Create models from saved options in model.info.'''
    Models = []
    for subdir in subdirs:
        print "Loading model from subdirectory: ", subdir
        Model = load_model(subdir,dryrun=dryrun)
        Models.append(Model)
    return Models

def new_models(subdirs,options):
    ''' Create new models with inputted options.'''
    Models = []
    for subdir in subdirs:
        Model = get_model(options)
        Models.append(Model)
    return Models
