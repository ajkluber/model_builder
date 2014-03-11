import os

import DMCModel
import HeterogeneousGoModel
import HomogeneousGoModel



'''
Author: Alexander Kluber

Purpose:
    A module that prepares the Model object which contains the functions and
parameters needed to prepare input files for coarse-grain protein simulations
in Gromacs. 

Description:
    models takes as input a dictionary of options or a subdirectory where to 
load options from a file. After checking options for consistency it returns 
it creates the corresponding Model objects.  

Changelog:
3-10-14 Created the options checker.
'''

def get_model_new(options):
    ''' Return a model with the following options. Options should have been
        previously checked for consistencty. '''
    type = options["Model_Code"]
    if type == "HomGo":
        model = HomogeneousGoModel.HomogeneousGoModel(disulfides=options["Disulfides"],
                                                    nonbond_param=options["nonbond_param"],
                                                    R_CD=options["R_CD"])
    elif type == "HetGo":
        model = HeterogeneousGoModel.HeterogeneousGoModel(disulfides=options["Disulfides"],
                                                    nonbond_param=options["nonbond_param"],
                                                    R_CD=options["R_CD"])
    elif type == "DMC":
        model = DMCModel.DMCModel()
    else:
        ## Due to previous error checking this should never be seen.
        print "ERROR! "
        print "Model: ",type," doesn't exist!"
    return model
    
def get_model(type):
    ''' SOON TO BE DEPRECATED. 3-10-14 AK'''
    if type == "HomGo":
        model = HomogeneousGoModel.HomogeneousGoModel()
    elif type == "HetGo":
        pass
        model = HeterogeneousGoModel.HeterogeneousGoModel()
    elif type == "DMC":
        model = DMCModel.DMCModel()
    else:
        print "ERROR. No such model exists."
    return model

def check_options(options):
    ''' Check that all options are compatible and in proper format. Any options
        that are omitted are given the value None (or 1 for nonbond_param) for
        completeness. Returns: options, a dictionary of completed options.'''

    ## List of supported models & corresponding representations.
    available_models = ["HomGo","HetGo","DMC"]
    beadmodels = {"HomGo":["CA"],"HetGo":["CA"]}
    contactopts = {"HetGo":["MJ","Bach","MC2004"]}

    modelcode = options["Model_Code"]
    beadmodel = options["Bead_Model"]
    
    print "Checking that options are consistent..."
    print "Inputted options:"
    print options, "\n"

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
    if options.has_key("Disulfides"):
        if options["Disulfides"] not in ["",None,False]:
            disulf = options["Disulfides"].split()
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

    ## Check for contact energies input.
    if options.has_key("Contact_Energies"):
        if options["Contact_Energies"] not in ["",None,False]:
            if contactopts.has_key(modelcode):
                if (options["Contact_Energies"] in contactopts[modelcode]):
                    contact_energies = options["Contact_Energies"]
                elif options["Contact_Energies"].split("/")[-1] == "BeadBead.dat":
                    if os.path.exists(options["Contact_Energies"]) == False:
                        print "ERROR!"
                        print "Contact energies option: ", options["Contact_Energies"], \
                            " points to a nonexistant file!"
                        print "Exiting."
                        raise SystemExit
                    else:
                        contact_energies = options["Contact_Energies"]
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
            contact_energies = None
    else:
        Contact_Energies = None
    options["Contact_Energies"] = contact_energies

    ## Check for R_CD option. This option fixes the ratio of contact (C) to
    ## dihedral (D) energy.
    if options.has_key("R_CD"):
        if options["R_CD"] not in ["",None,False]:
            try:
                R_CD = float(options["R_CD"])
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
    if options.has_key("nonbond_param"):
        if options["nonbond_param"] not in ["",None,False]:
            try:
                nonbond_param = float(options["nonbond_param"])
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
    ## using heavy atom contacts within a cutoff of a given residue.
    if options.has_key("Cutoff"):
        if options["Cutoff"] not in ["",None,False]:
            try:
                cutoff = float(options["Cutoff"])
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
    if options.has_key("Solvent"):
        if options["Solvent"] in ["",None,False]:
            solvent = None
        else: 
            print "Error! Solvent option not implemented!"
            print "Solvent variable: ", options["Solvent"], " not allowed. Only: ",["",None]
            print "Exiting."
            raise SystemExit
    else:
        solvent = None
    options["Solvent"] = solvent


    print "Options cleared!"
    print "Using options:"
    print options
    print "Proceeding...\n"
            
    return options

def load_model(path):
    ''' Given path that contains model.info options file. Read in options and
        create corresponding model.'''
    info_file = open(path+'model.info','r')
    line = info_file.readline()
    while line != '':
        field = line.split()[1]
        value = info_file.readline()
        options[field] = value[:-1]
    options = check_options(options)
    Model = get_model_new(options)
    return Model

def load_models(subdirs):
    ''' Create models from saved options in model.info.'''
    Models = []
    for subdir in subdirs:
        Model = load_model(subdir)
        Models.append(Model)
    return Models

def new_models(subdirs,options):
    ''' Create new models with options.'''
    Models = []
    for subdir in subdirs:
        Model = get_model_new(options)
        Models.append(Model)
    return Models
