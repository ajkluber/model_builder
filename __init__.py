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

    modelcode = options["Model_Code"]
    beadmodel = options["Bead_Model"]
    
    print "Checking that options are consistent..."
    print "Inputted options:"
    print options, "\n"

    if modelcode not in available_models:
        print "ERROR! Invalid Model_Code"
        print "Model: ", modelcode, "  not within available models: ",available_models
        print "Exiting."
        raise SystemExit

    if beadmodel not in beadmodels[modelcode]:
        print "ERROR! Invalid Bead_Model"
        print "Bead model: ", beadmodel, " not available for model: ", modelcode
        print "Model: ",modelcode," has the following bead models: ", beadmodels[modelcode]
        print "Exiting."
        raise SystemExit
    
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
