''' Check inputs to SmogCalpha '''

import numpy as np
import os
import shutil

import models.SmogCalpha as SmogCalpha


global negvals
negvals = ["None",None,"",False,"False"]

def get_pairwise_params(pairwise_params_file,model_params_file):
    ''' Grab pairwise_params from file. '''
    model_param_values = np.loadtxt(model_params_file)

    p_lines = [ x.rstrip("\n") for x in open(pairwise_params_file,"r").readlines() ]

    pairs = []
    pairwise_param_assignment = []
    pairwise_type = [] 
    pairwise_other_params = []

    for p in p_lines[1:]:
        data = p.split() 
        pairs.append([int(data[0]),int(data[1])])
        pairwise_param_assignment.append(int(data[2]))
        pairwise_type.append(int(data[3]))
        temp = []
        for otherparam in data[4:]:
            temp.append(float(otherparam))
        pairwise_other_params.append(temp)

    pairs = np.array(pairs) 
    pairwise_param_assignment = np.array(pairwise_param_assignment)
    pairwise_type = np.array(pairwise_type)

    return pairs,pairwise_param_assignment,model_param_values,pairwise_type,pairwise_other_params

def check_contact_args(inputs,pairs_file,pairwise_params_file,model_params_file,epsilonbar):
    ''' Check input arguments for pairwise interactions'''
    pairs = None
    epsilon_bar = None
    inputs["Defaults"] = True

    ## Load pairs or pairs+parameters if given.
    if pairwise_params_file in negvals:
        if pairs_file in negvals:
            print "ERROR! specify either --pairs <filename>  or --pairwise_params_file <filename>!"
            print "Exiting"
            raise SystemExit
        else:
            if not os.path.exists(pairs_file):
                print "ERROR! %s doesn't exist! " % pairs_file
                print "Exiting"
                raise SystemExit
            else:
                pairs = np.loadtxt("%s" % pairs_file,dtype=int)
        ##For removing the unnecessary columns from the smog contact map output
        if len(pairs[0,:]) == 4:
            pairs = pairs[:,np.array([1,3])]
    else:
        pairs,pairwise_param_assignment,model_param_values,pairwise_type,pairwise_other_params = get_pairwise_params(pairwise_params_file,model_params_file)
        inputs["pairwise_param_assignment"] = pairwise_param_assignment
        inputs["model_param_values"] = model_param_values 
        inputs["pairwise_type"] = pairwise_type
        inputs["pairwise_other_parameters"] = pairwise_other_params
        inputs["pairwise_params_file_location"] = pairwise_params_file
        inputs["model_params_file_location"] = model_params_file
        inputs["Defaults"] = False

    inputs["Contacts"] = pairs
    ## Check if average contact strength parameter, epsilon_bar, is set. This 
    ## keeps the average contact strength normalized to some number to maintain
    ## the same stability if the parameters are modified.
    if epsilonbar not in negvals:
        try:
            epsilon_bar = float(epsilonbar)
        except:
            print "TypeError! epsilon_bar value must be a float!"
            print "Exiting."
            raise SystemExit
        inputs["Epsilon_Bar"] = epsilon_bar

    return inputs

def check_fitting_args(inputs,fittingdata,fittingincludes,fittingsolver,fittingallowswitch,fittingparamsfile):
    ''' Check parameter fitting options.

    If parameter fitting is being used check that the fitting inputs make sense:
      fitting_data          indicates the type of data to fit
      fitting_includes      allows fitting over multiple subdirectories.
      fitting_solver        choses the algorithm to select a solution.
      fitting_allowswitch   specifies if interactions are allowed to change from attractive/repulsive.
      fitting_params        specifies which parameters are being fit.
    '''
    fittingdatas = ["ddG_MC2004","RMSF","FRET","contact_Qi"]
    fittingallowswitches = ["True","False"]
    fittingsolvers = ["Levenberg","TSVD","TSVD_Cplex"]

    if fittingdata in negvals:
        fitting_data = None
        fitting_includes = [ None ]
    else:
        if fittingdata in fittingdatas:
            fitting_data = fittingdata
        else:
            print "KeyError! Fitting_Data must be one of:", fittingdatas
            print "Exiting."
            raise SystemExit
        if fittingincludes in negvals:
            fitting_includes = [ inputs["PDB"].split(".pdb")[0] ]
        else:
            fitting_includes = [ fittingincludes[i].split(".pdb")[0] for i in range(len(fittingincludes)) ]

    if fittingallowswitch in fittingallowswitches:
       fitting_allowswitch = fittingallowswitch 
    else:
        fitting_allowswitch = "False"

    if fittingsolver not in fittingsolvers:
        print "ERROR! Fitting_Solver must be one of:", fittingsolvers
        print "Exiting."
        raise SystemExit
    else:
        fitting_solver = fittingsolver

    if fittingparamsfile in negvals:
        fitting_params = None
        fitting_params_file = None
    else:
        if not os.path.exists(fittingparamsfile):
            print "ERROR! Fitting_Params file doesn't exist:", fittingparamsfile
            print "Exiting."
            raise SystemExit
        else:
            fitting_params = np.loadtxt(fittingparamsfile,dtype=int)
            fitting_params_file = fittingparamsfile

    inputs["Fitting_Data"] = fitting_data
    inputs["Fitting_Includes"] = fitting_includes
    inputs["Fitting_Solver"] = fitting_solver
    inputs["Fitting_AllowSwitch"] = fitting_allowswitch
    inputs["Fitting_Params"] = fitting_params
    inputs["Fitting_Params_File"] = fitting_params_file

    return inputs

def check_disulfide_args(inputs,inputdisulfides):
    ''' Check disulfide argument '''
    if inputdisulfides in negvals:
        disulfides = None
    else:
        disulf = inputdisulfides
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
    inputs["Disulfides"] = disulfides
    return inputs

def new_args(args):
    ''' Check new input arguments '''
    available_models = ["HomGo","HetGo","DMC"]
    beadmodels = {"HomGo":["CA"],"HetGo":["CA"]}
    fittingdatas = ["ddG_MC2004","RMSF","FRET","contact_Qi"]
    inputs = {}

    print "Checking input options..."
    ## Check for pdb existence
    for pdb in args.pdbs:
        if not os.path.exists(pdb):
            print "ERROR! pdb %s doesn't exist! " % pdb
            print "Exiting"
            raise SystemExit
    model_code = args.model_code
    beadmodel = args.bead_model 

    inputs["Iteration"] = 0
    inputs["PDB"] = args.pdbs[0]
    inputs["Model_Code"] = args.model_code
    inputs["Bead_Model"] = args.bead_model
    ## Check if model code is legal.
    if model_code not in available_models:
        print "ERROR! Invalid Model_Code"
        print "Model: ", model_code, "  not within available models: ",available_models
        print "Exiting."
        raise SystemExit

    ## Check if bead model is available for model code.
    if beadmodel not in beadmodels[model_code]:
        print "ERROR! Invalid Bead_Model"
        print "Bead model: ", beadmodel, " not available for model: ", model_code
        print "Model: ",model_code," has the following bead models: ", beadmodels[model_code]
        print "Exiting."
        raise SystemExit

    ## Check all contact-related inputs
    inputs = check_contact_args(inputs,args.pairs,args.pairwise_params_file,args.model_params_file,args.epsilon_bar)

    ## Check parameter fitting inputs
    inputs = check_fitting_args(inputs,args.fitting_data,args.fitting_includes,args.fitting_solver,args.fitting_allowswitch,args.fitting_params_file)

    ## Check disulfide list
    inputs = check_disulfide_args(inputs,args.disulfides)

    ## Dry run flag will prevent any simulations from being submitted. Used to
    ## see if file preparation runs smoothly.
    if args.dry_run == True:
        dryflag = True
    else:
        dryflag = False
    inputs["Dry_Run"] = dryflag

    print "Using model options:"
    keys = inputs.keys()
    keys.sort()
    for key in keys:
        if key in ["Contacts","model_param_values","pairwise_other_parameters",
                   "pairwise_type","pairwise_param_assignment","fitting_params",
                    "Fitting_Params"]:
            if inputs[key] == None:
                print "  ", key , " = ", inputs[key]
            else:
                print "  ", key , " = [", inputs[key][0], "...",inputs[key][-1], "]"
        else:
            print "  ", key , " = ", inputs[key]
    return inputs


def load_args(subdir,dry_run):

    info_file = open('%s/model.info' % subdir,'r')
    line = info_file.readline()
    inputs = {"Dry_Run":dry_run}
    while line != '':
        field = line.split()[1]
        value = info_file.readline()
        if field == "Reference":
            break
        elif field in ["Interaction_Groups","Model_Name",
                        "Backbone_params","Backbone_param_vals",
                        "Interaction_Types","Tf_iteration","Mut_iteration",
                        "Contact_Type","Contact_params","Contact_Energies"]:
            pass
        elif field == "Iteration":
            inputs[field] = int(value.rstrip("\n"))
        elif field == "N_Native_Contacts":
            if value.rstrip("\n") != "None":
                inputs[field] = int(value.rstrip("\n"))
            else:
                inputs[field] = None
        else:
            inputs[field] = value.rstrip("\n")
        line = info_file.readline()

    if not os.path.exists("%s/pairs.dat" % subdir):
        shutil.copy("%s/contacts.dat" % subdir,"%s/pairs.dat" % subdir)

    pairs_file = "%s/pairs.dat" % subdir

    pairwise_params_file = inputs["Pairwise_Params_File"]
    model_params_file = inputs["Model_Params_File"]
    epsilonbar = inputs["Epsilon_Bar"] 
    fittingdata = inputs["Fitting_Data"]
    fittingincludes = inputs["Fitting_Includes"].split()
    fittingallowswitch = inputs["Fitting_AllowSwitch"]
    fittingsolver = inputs["Fitting_Solver"]
    fittingparamsfile = inputs["Fitting_Params_File"]
    disulfides = inputs["Disulfides"]
        
    ## Check all contact-related inputs
    inputs = check_contact_args(inputs,pairs_file,pairwise_params_file,model_params_file,epsilonbar)

    ## Check parameter fitting inputs
    inputs = check_fitting_args(inputs,fittingdata,fittingincludes,fittingsolver,fittingallowswitch,fittingparamsfile)

    ## Check disulfide list
    inputs = check_disulfide_args(inputs,disulfides)

    print "Using model options:"
    keys = inputs.keys()
    keys.sort()
    for key in keys:
        if key in ["Contacts","model_param_values","pairwise_other_parameters",
                   "pairwise_type","pairwise_param_assignment","fitting_params",
                    "Fitting_Params"]:
            if inputs[key] == None:
                print "  ", key , " = ", inputs[key]
            else:
                print "  ", key , " = [", inputs[key][0], "...",inputs[key][-1], "]"
        else:
            print "  ", key , " = ", inputs[key]
    return inputs

def load_model(subdir,dry_run=False):
    ''' Read model.info files in subdirectories and create models.'''

    options = load_args(subdir,dry_run)
    model = SmogCalpha.SmogCalpha(**options)

    return model

def load_models(subdirs,dry_run=False):
    ''' Create models from saved options in model.info.'''
    Models = []
    for subdir in subdirs:
        print "Loading model from subdirectory: ", subdir
        #Model = load_model(subdir,dry_run=dry_run)
        options = load_args(subdir,dry_run)
        Model = SmogCalpha.SmogCalpha(**options)
        Models.append(Model)
    return Models

def new_models(subdirs,options):
    ''' Create new models with inputted options.'''
    Models = []
    for subdir in subdirs:
        options["PDB"] = subdir+".pdb"
        model = SmogCalpha.SmogCalpha(**options)

        Models.append(model)
    return Models
