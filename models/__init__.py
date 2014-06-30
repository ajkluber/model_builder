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

import numpy as np
import os

import DMCModel
import HeterogeneousGoModel
import HomogeneousGoModel
import SecondaryTertiaryGoModel
import CalphaBase
import SmogCalpha

def get_contact_params(paramfile):
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

def check_options(inputoptions,firstpass=False):
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
                    contacts = None
                    contact_epsilons = None
                    contact_deltas = None
                elif (inputoptions["Contact_Energies"].endswith(".dat")) or (inputoptions["Contact_Energies"].endswith(".params")):
                    if not os.path.exists(inputoptions["Contact_Energies"]):
                        print "ERROR!"
                        print "Contact energies option: ", inputoptions["Contact_Energies"], \
                            " points to a nonexistant file!"
                        print "Exiting."
                        raise SystemExit
                    else:
                        contacts,contact_epsilons,contact_deltas = get_contact_params(inputoptions["Contact_Energies"])
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
    options["Contact_Energies"] = contact_energies
    options["Contacts"] = contacts
    options["Contact_Epsilons"] = contact_epsilons
    options["Contact_Deltas"] = contact_deltas

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

    ## Dry run flag will prevent any simulations from being submitted. Used to
    ## see if file preparation runs smoothly.
    if inputoptions["Dry_Run"] == True:
        dryflag = True
    else:
        dryflag = False
    options["Dry_Run"] = dryflag

    if firstpass == False:
        print "Model options cleared! Using model options:"
        for key in options.keys():
            print "  ", key , " = ", options[key]
            
    return options

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
            disulfides=options["Disulfides"],
            modelcode=options["Model_Code"],
            contact_energies=options["Contact_Energies"],
            Tf_iteration=options["Tf_Iteration"],
            Mut_iteration=options["Mut_Iteration"],
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
                disulfides=options["Disulfides"],
                modelcode=options["Model_Code"],
                contact_energies=options["Contact_Energies"],
                Tf_iteration=options["Tf_Iteration"],
                Mut_iteration=options["Mut_Iteration"],
                dryrun=options["Dry_Run"])

        Models.append(model)
    return Models
