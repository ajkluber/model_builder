""" Check inputs to 


"""

import numpy as np
import os

import models.SmogCalpha as SmogCalpha

def get_contact_params(paramfile,contact_type):
    """ Read contact info from inputted contact parameter file """
    if paramfile.endswith(".dat"):
        beadbead = np.loadtxt(paramfile,dtype=str) 
        contacts = []
        contact_epsilons = []
        contact_deltas = []
        for i in range(len(beadbead[:,4])):
            if beadbead[i,4] not in ["0","ss"]:
                contacts.append(beadbead[i,0:2].astype(int))
                contact_epsilons.append(beadbead[i,6].astype(float))
                contact_deltas.append(beadbead[i,7].astype(float))
        contacts = np.array(contacts)
        contact_epsilons = np.array(contact_epsilons)
        contact_deltas = np.array(contact_deltas)
    elif paramfile.endswith(".params"):
        params = np.loadtxt(paramfile)
        contacts = params[:,0:2].astype(int)
        contact_epsilons = params[:,2].astype(float)
        contact_deltas = params[:,3].astype(float)
    else:
        print "ERROR!"
        print " Unsupported file type:", paramfile
        print " Provide a .dat or .params file"
        print " Exiting."
        raise SystemExit

    return contacts,contact_epsilons,contact_deltas

def new_args(args):
    """ Check new input arguments """

    contacttypes = ["LJ1210","Guassian"]
    inputs = {}
    ## To Do:
    ## Check for
    ##  - existence of pdbs file 
    ##  - check existence of contacts file, load contacts
    ##  - 
    ## 
    ## 

    for pdb in args.pdbs:
        if not os.path.exists(pdb):
            print "ERROR! pdb %s doesn't exist! " % pdb
            print "Exiting"
            raise SystemExit

    inputs["PDB"] = args.pdbs[0]

    if args.contact_type == "None":
        contact_type = "LJ1210"
    else:
        if args.contact_type not in contacttypes:
            print "ERROR! %s not recognized. --contact_type must be from: " % args.contacts, contacttypes
            print "Exiting"
            raise SystemExit
        else:
            contact_type = args.contact_type

    if (args.contacts == "None") and (args.contact_params == "None"):
        print "ERROR! specify either --contacts <filename>  or --contact_params <filename>!"
        print "Exiting"
        raise SystemExit
    elif (args.contacts == "None") and (args.contact_params != "None"):
        get_contact_params(paramfile,contact_type)
    elif (args.contacts != "None") and (args.contact_params == "None"):
        if not os.path.exists(args.contacts):
            print "ERROR! %s doesn't exist! " % args.contacts
            print "Exiting"
            raise SystemExit
        else:
            contacts = np.loadtxt("%s" % args.contacts,dtype=int)


    pass

def check_new_contact_inputs(args):
    pass
    #if args.contacts == "None":

def load_model_info():
    pass

def check_options(inputoptions,firstpass=False):
    ''' Check that all options are compatible and in proper format. Any options
        that are omitted are given the value None (or 1 for nonbond_param) for
        completeness. Returns: options, a dictionary of completed options.'''

    ## List of supported models & corresponding representations.
    available_models = ["HomGo","HetGo","DMC"]
    beadmodels = {"HomGo":["CA"],"HetGo":["CA"]}
    contactopts = {"HetGo":["MJ","Bach","MC2004","FRETFit","RMSFit","SecTer"]}
    contacttypes = ["LJ1210","Guassian"]
    fittingopts = ["ddG_MC2004","RMSF","FRET"]

    modelcode = inputoptions["Model_Code"]
    beadmodel = inputoptions["Bead_Model"]
    
    print "Checking that model options are consistent..."
    #print "Inputted model options:"
    #for key in inputoptions.keys():
    #    print "  ", key , " = ", inputoptions[key]

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

    ## Check disulfide options.
    disulfides = check_disulfide_options(inputoptions)
    options["Disulfides"] = disulfides

    ## Check contact type and epsilon_bar
    contact_type, epsilon_bar = check_contact_type_and_epsilon_bar(inputoptions,contacttypes)
    options["Contact_Type"] = contact_type
    options["Epsilon_Bar"] = epsilon_bar

    ## Check options for contacts
    contact_energies, contacts, contact_epsilons, contact_deltas = check_contact_params_options(inputoptions,modelcode,contactopts,contact_type)
    options["Contact_Energies"] = contact_energies
    options["Contacts"] = contacts
    options["Contact_Epsilons"] = contact_epsilons
    options["Contact_Deltas"] = contact_deltas


    ## Get contacts if not already specified
    print inputoptions
    if (contacts == None):
        if firstpass:
            contacts = None
        else:
            if inputoptions.has_key("Contacts"):
                if inputoptions["Contacts"] in ["None",None]:
                    contacts = None
                else:
                    if not os.path.exists(inputoptions["Contacts"]):
                        print "ERROR!"
                        print "Contact energies option: ", inputoptions["Contacts"], \
                            " points to a nonexistant file!"
                        print "Exiting."
                        raise SystemExit
                    else:
                        contacts = np.loadtxt(inputoptions["Contacts"],dtype=int)
            else:
                print "Error! Contacts option must be given if contact params not set with contact_energies!"
                print "Exiting."
                raise SystemExit
        options["Contacts"] = contacts

    ## Check if procedural indicator Tf_Iteration is set
    if inputoptions.has_key("Tf_Iteration"):
        try:
            Tf_iteration = int(inputoptions["Tf_Iteration"])
        except:
            print "TypeError! Tf_iteration value must be a int!"
            print "Exiting."
            raise SystemExit
    else:
        Tf_iteration = 0
    options["Tf_Iteration"] = Tf_iteration

    ## Check if procedural indicator Mut_Iteration is set
    if inputoptions.has_key("Mut_Iteration"):
        try:
            Mut_iteration = int(inputoptions["Mut_Iteration"])
        except:
            print "TypeError! Mut_iteration value must be a int!"
            print "Exiting."
            raise SystemExit
    else:
        Mut_iteration = 0
    options["Mut_Iteration"] = Mut_iteration

    ## Fitting_Data indicates the type of data to fit
    if inputoptions.has_key("Fitting_Data"):
        if inputoptions["Fitting_Data"] in fittingopts:
            Fitting_Data = inputoptions["Fitting_Data"]
        else:
            print "KeyError! Fitting_Data must be one of:", fittingopts
            print "Exiting."
            raise SystemExit
    else:
        Fitting_Data = None
    options["Fitting_Data"] = Fitting_Data

    ## Fitting_Includes allows fitting over multiple subdirectories. Defaults
    ## to only the pdb given.
    if inputoptions.has_key("Fitting_Includes"):
        if len(inputoptions["Fitting_Includes"]) > 1:
            Fitting_Includes = inputoptions["Fitting_Includes"]
        else:
            Fitting_Includes = [ inputoptions["PDB"].split(".pdb")[0] ]
    else:
        Fitting_Includes = [ inputoptions["PDB"].split(".pdb")[0] ]
    options["Fitting_Includes"] = Fitting_Includes

    ## Dry run flag will prevent any simulations from being submitted. Used to
    ## see if file preparation runs smoothly.
    if inputoptions["Dry_Run"] == True:
        dryflag = True
    else:
        dryflag = False
    options["Dry_Run"] = dryflag

    if not firstpass:
        print "Model options cleared! Using model options:"
        for key in options.keys():
            if key in ["Contacts","Contact_Epsilons","Contact_Deltas"]:
                if options[key] == None:
                    print "  ", key , " = ", options[key]
                else:
                    print "  ", key , " = ", options[key][:3], "...",options[key][-3:]
            else:
                print "  ", key , " = ", options[key]
            
    return options

def check_disulfide_options(inputoptions):
    """ Check that disulfides option is consistent 

    Description:

        Check for disulfides. If yes disulfids should be in format (e.g.): 
    [ 1, 6, 7, 20] where the first pair of numbers 1,6 correspond to a 
    disulfide bond between those residue indices according to clean.pdb indices.
    """

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
    return disulfides

def check_contact_params_options(inputoptions,modelcode,contactopts,contact_type):
    """ Check the contact params options for consistency 

    Description:

        Check for contact energies input. Only valid for certain models. Format
    for contact_energies should be an allowed option for the model (e.g. 'MJ'
    for 'HetGo') or a path to a Beadbead.dat with previously saved contact
    energies.
    """

    
    if inputoptions.has_key("Contact_Energies"):
        if inputoptions["Contact_Energies"] not in ["","None",None,False]:
            if contactopts.has_key(modelcode):
                if (inputoptions["Contact_Energies"] in contactopts[modelcode]):
                    contact_energies = inputoptions["Contact_Energies"]
                    contacts = None
                    contact_epsilons = None
                    contact_deltas = None
                elif (inputoptions["Contact_Energies"].endswith(".dat")) \
                    or (inputoptions["Contact_Energies"].endswith(".params")) \
                    or (inputoptions["Contact_Energies"].endswith(".contacts")):
                    if not os.path.exists(inputoptions["Contact_Energies"]):
                        print "ERROR!"
                        print "Contact energies option: ", inputoptions["Contact_Energies"], \
                            " points to a nonexistant file!"
                        print "Exiting."
                        raise SystemExit
                    else:
                        contacts,contact_epsilons,contact_deltas = get_contact_params(inputoptions["Contact_Energies"],contact_type)
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
                contacts = None
                contact_epsilons = None
                contact_deltas = None
    else:
        contact_energies = None
        contacts = None
        contact_epsilons = None
        contact_deltas = None
    return contact_energies, contacts, contact_epsilons, contact_deltas

def check_contact_type_and_epsilon_bar(inputoptions,contacttypes):
    """ Check contact type option and epsilon bar option """

    if inputoptions.has_key("Contact_Type"):
        if inputoptions["Contact_Type"] not in ["","None",None,False]:
            if inputoptions["Contact_Type"] not in contacttypes:
                print "KeyError! contact_type be in ", contacttypes
                print "Exiting."
                raise SystemExit
            else:
                contact_type = inputoptions["Contact_Type"]
        else:
            contact_type = "LJ1210"
    else:
        contact_type = "LJ1210"

    ## Check if average contact strength parameter, epsilon_bar, is set. This 
    ## keeps the average contact strength normalized to some number to maintain
    ## the same stability if the parameters are modified.
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
    return contact_type, epsilon_bar


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
    options["PDB"] = subdir+".pdb"

    model = SmogCalpha.SmogCalpha(
            options["PDB"],
            contacts=options["Contacts"],
            contact_epsilons=options["Contact_Epsilons"],
            contact_deltas=options["Contact_Deltas"],
            epsilon_bar=options["Epsilon_Bar"],
            contact_type=options["Contact_Type"],
            disulfides=options["Disulfides"],
            modelcode=options["Model_Code"],
            contact_energies=options["Contact_Energies"],
            Tf_iteration=options["Tf_Iteration"],
            Mut_iteration=options["Mut_Iteration"],
            fitting_data=options["Fitting_Data"],
            fitting_includes=options["Fitting_Includes"],
            dryrun=options["Dry_Run"])

    return model

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
        options["PDB"] = subdir+".pdb"
        model = SmogCalpha.SmogCalpha(
                options["PDB"],
                contacts=options["Contacts"],
                contact_epsilons=options["Contact_Epsilons"],
                contact_deltas=options["Contact_Deltas"],
                epsilon_bar=options["Epsilon_Bar"],
                contact_type=options["Contact_Type"],
                disulfides=options["Disulfides"],
                modelcode=options["Model_Code"],
                contact_energies=options["Contact_Energies"],
                Tf_iteration=options["Tf_Iteration"],
                Mut_iteration=options["Mut_Iteration"],
                fitting_data=options["Fitting_Data"],
                fitting_includes=options["Fitting_Includes"],
                dryrun=options["Dry_Run"])

        Models.append(model)
    return Models
