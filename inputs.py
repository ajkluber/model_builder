""" Check inputs for making a CoarseGrainedModel """

import numpy as np
import re
import os
import shutil
import ConfigParser
import mdtraj

from models import StructureBasedModel as SBM

#############################################################################
# Helper functions to load in models from .ini files
#############################################################################

def load_model(name,dry_run=False):
    """ Load a model from a config file.
    
    Constructs a models.Model class or one of its subclasses using the 
    provided config file. 
    
    All files handled by load_model is assumed to be formatted with 
    their respective outside use formats. I.E. Residues numbering starts 
    from 1, not zero.
    All inputs to methods should be assumed to take the python way of 
    numbering, i.e. start from zero. 
    By Design, this inputs.py file should be where any and all 
    conversions take place for consistency. The assumed conversions and 
    defaults are listed under Default Assumptions. Some conversions 
    might still be necessary inside other methods, but should generally 
    be avoided.
    
    Args:
        name(string): Full name of the config file to load
        
    Return:
        Model: A Model object constructed from the config file
        Dict.: List of fitting options
        
    Default Assumptions:
        pairs start from 1 (residue-residue contacts for CA-CA 
            potentials). Get converted to starting from zero when 
            calling model.add_pairs
        Default_Energy_parameters for angle bonds are assumed to be 
            given in the gromacs way (/rad^2). Inside the model_builder 
            package, the energy and angles are also assumed to use 
            radians. They are converted accordingly in the outputted 
            files. Therefore, all conversions should only take place in 
            the outputs package.
    
    """
    
    #Load options and construct topology
    modelopts, fittingopts = load_config(name)
    modelopts["dry_run"] = dry_run
    top = mdtraj.load(modelopts["topology"])
    model = SBM(top.topology, bead_repr=modelopts["bead_repr"])
    
    #load a reference set to base a model off of
    if modelopts["reference"] is None:
        traj = top
    else:
        traj = mdtraj.load(modelopts["reference"])
    model.set_reference(traj)
    
    ##add backbone and disulfides
    model.assign_backbone()
    #add disulfides
    if modelopts["disulfides"] is not None:
        disulf = modelopts["disulfides"]
        disulfides = []
        for i in range(len(disulf)/2):
            disulfides.append([disulf[2*i]-1, disulf[2*i+1]-1])
        model.assign_disulfides(disulfides)
    model.add_sbm_backbone()
    
    #Check for pair options
    if modelopts["pairs"] is None:
        if modelopts["pairwise_params_file"] is None:
            model.assign_contacts()
            model.add_sbm_contacts
        else:
            #use the modelopts["pairwise_params"] file
            pairs, pairs_index_number, pairs_potential_type, pairs_args = parse_pairwise_params(modelopts["pairwise_params_file"])
            model.add_pairs(pairs)
            pairopts = []
            if modelopts["model_params_file"] is None:
                print "Warning: No model_params_file specified. Defaulting all epsilons to 1"
                epsilons = np.ones(np.shape(pairs)[0])
            else:
                if os.path.exists(modelopts["model_params_file"]):
                    epsilons = np.loadtxt(modelopts["model_params_file"], comments="#")
                else:
                    raise IOError("{} does not exist!".format(modelopts["model_params_file"]))

            if not np.shape(epsilons)[0] == len(pairs):
                raise IOError("Number of model params not equal to number of pairwise params")
                
            for i in range(len(pairs)):
                code = pairs_potential_type[i]
                atm1, atm2 = model.mapping._contact_pairs[i]
                eps = epsilons[i]
                pot_type = pairs_potential_type[i]
                opts = [pot_type, atm1, atm2, eps]
                for args in pairs_args[i]:
                    opts.append(args)
                pairopts.append(opts)
            model.Hamiltonian._add_pairs(pairopts)
    else:
        #use the modelopts["pairs"] file
        model.add_pairs(np.loadtxt(modelopts["pairs"]).astype(int) - 1)
        model.add_sbm_contacts()
        
    #Assign epsilons for fitting. Default: All pair interactions
    if modelopts["parameters_to_fit_file"] is None:
        parameters_to_fit = np.arange(np.shape(model.Hamiltonian._pairs)[0])
    else:
        parameters_to_fit = np.loadtxt(modelopts["parameters_to_fit_file"], comments="#").astype(int)
    
    model.assign_fitted_epsilons(parameters_to_fit)
    
    return model,fittingopts

def load_models(names,dry_run=False):
    """Create models from saved options in <name>.ini"""
    Models = []
    Fittingopts = []
    for name in names:
        print "Loading model for: %s" % name
        model,fittingopts = load_model(name,dry_run=dry_run)
        Models.append(model)
        Fittingopts.append(fittingopts)
    return Models,Fittingopts

def load_config(name):
    """Parse options from <name>.ini file"""
     
    if not os.path.exists("%s" % name):
        raise IOError("%s doesn't exist!" % name)
    config = ConfigParser.SafeConfigParser(allow_no_value=True)
    config.read("%s" % name)

    # Populate all fields as None 
    modelopts = _empty_model_opts()
    fittingopts = _empty_fitting_opts()

    print "Creating model according to %s.ini" % name
    print "Options not shown default to None"
    load_model_section(config.items("model"),modelopts)
    load_fitting_section(config,modelopts,fittingopts)
    #_add_pair_opts(modelopts) 
    return modelopts,fittingopts

def load_model_section(modelitems,modelopts):
    """Parse [model] options from config .ini file"""
    
    #acceptable strings for boolean
    bool_valid_check = ["true", "True", "1", "T", "t"] 
    
    #generic checks/handling of certain options
    check_boolean = ["using_sbm_gmx", "umbrella", "simple_disulfides"]
    check_exists = ["topology", "starting_gro", "reference"]
    
    print "Model options:"
    for item,value in modelitems:
        if value in [None,""]:
            pass
        else:
            print "  %-20s = %s" % (item,value)
            #Generic Checks
            if item in check_boolean:
                value = value in bool_valid_check
            if item in check_exists:
                if not os.path.exists("%s" % value):
                    raise IOError("%s file does not exist! Check config file inputs" % value)
            
            #specific checks
            if item == "n_native_pairs":
                value = int(value)
            elif item == "epsilon_bar":
                value = float(value)
            
            elif item == "disulfides":
                value = [ int(x) for x in re.split(",\s+|\s+", value.strip("[ | ]"))]
                if (len(value) % 2) != 0:
                    raise IOError("len(disulfides) should be even. Invalid input: %s " % value.__repr__())
            elif item.endswith("_file"):
                if not os.path.exists(value):
                    raise IOError("%s file does not exist! Check config file inputs" % value)
            elif item == "cb_volume":
                if value.endswith(".dat"):
                    if not os.path.exists(value):
                        raise IOError("%s file does not exist! Check config file inputs" % value)
                else:
                    if value not in ["average","flavored"]:
                        raise IOError("cb_volume value must be: average, flavored, or filename")
            elif item == "backbone_param_vals":
                value = eval(value)
            modelopts[item] = value

    if modelopts["bead_repr"] == "CACB" and (modelopts["cb_volume"] not in ["average","flavored"]):
        raise IOError("If bead_repr = CACB, then must set cb_volume to: average or flavored")

def load_fitting_section(config,modelopts,fittingopts):
    """Parse [fitting] options from config .ini file"""
    # special fitting checks is for package specific options
    # assigns based on keys, functions should be at end of file
    special_fitting_checks = {"FRET":FRET_fitopts_load, "tmatrix":tmatrix_fitopts_load}
    bool_valid_check = ["true", "True", "1", "T", "t"]
    
    if config.has_section("fitting"):
        if config.has_option("fitting","data_type") and (config.get("fitting","data_type") in special_fitting_checks):
            checkfunction = special_fitting_checks[config.get("fitting","data_type")]
            check_special = True
        else:
            check_special = False

        print "\nFitting options:"
        for item,value in config.items("fitting"):
            if value in [None,""]:
                pass
            else:
                print "  %-20s = %s" % (item,value)
                if item == "iteration":
                    value = int(value)
                elif item == "include_dirs":
                    value = value.split()
                elif item == "allow_switch":
                    value = value in bool_valid_check
                elif item == "constrain_avg_eps":
                    value = value in bool_valid_check
                elif item == "nonnative":
                    value = value in bool_valid_check
                elif item in ["equil_walltime","walltime"]:
                    if re.match("\d\d:\d\d:\d\d",value) == None:
                        raise IOError(" %s must be a time in format HH:MM:SS" % item)
                elif item == "parameters_to_fit":
                    if not os.path.exists(value):
                        raise IOError("%s file does not exist! Check config file inputs" % value)
                    else:
                        modelopts["fitting_params"] = np.loadtxt(value,dtype=int)
                elif item == "cutoffs":
                    value = [ float(x) for x in re.split(",\s+|\s+", value.strip("[ | ]"))] 
                elif item == "simplify_lambdas":
                    value = value in bool_valid_check
                elif check_special:
                    value = checkfunction(item, value)
                fittingopts[item] = value
                        
#############################################################################
# Internal functions to load in models from .ini files
#############################################################################
def _empty_fitting_opts():
    """Fitting options to check for"""
    opts = ["data_type","include_dirs","solver",
            "iteration","n_processors","equil_walltime",
            "allow_switch","parameters_to_fit",
            "nonnative","last_completed_task","cutoffs","simplify_lambdas"]         
    fittingopts = { opt:None for opt in opts }
    fittingopts["simplify_lambdas"] = False
    return fittingopts

def _empty_model_opts():
    """Model options to check for"""
    opts = ["pairs_file","pairwise_params_file",
            "model_params_file","epsilon_bar","starting_gro",
            "defaults","bead_repr","cb_volume","disulfides",
            "simple_disulfides","n_native_pairs","contact_type",
            "pairs","pairwise_other_parameters",
            "pairwise_param_assignment","n_processors",
            "pairwise_type","verbose","dry_run",
            "using_sbm_gmx","umbrella", "topology", "reference", "parameters_to_fit_file"] 
    modelopts = { opt:None for opt in opts }
    return modelopts

def parse_pairwise_params(pairwise_file):
    """ parse the pairwise_params file and output necessary values"""
    
    fopen = open(pairwise_file, "r")
    pairs = []
    pairs_index_number = []
    pairs_potential_type = []
    pairs_args = []
    key={"8":"LJ12GAUSSIAN", "4": "GAUSSIAN", "2":"LJ1210", "5":"LJ12TANHREP"}
    wordkeys = [key[i] for i in key]
    wordkeys.append("LJ12GAUSSIANTANH")
    count = 0
    for line in fopen:
        count += 1
        data = line.strip().split()
        if not data[0][0] is "#":
            #convert to pythonic indices
            pairs.append([int(data[0])-1, int(data[1])-1]) 
            pairs_index_number.append(int(data[2]))
            
            if data[3] in key:
                #using number keys for functions, convert to Word key
                pot_type = key[data[3]] 
            elif data[3] in wordkeys:
                #check if using word keys
                pot_type = data[3]
            else:
                print "Unknown function type for param: %s %s\n" %(data[0], data[1])
                print "See line %d"
            
            pairs_potential_type.append(pot_type)
            
            pairs_args.append([float(val) for val in data[4:]])
    
    return pairs, pairs_index_number, pairs_potential_type, pairs_args
    
def _add_pairwise_params(modelopts):
    """Parse pairwise_params file"""

    model_param_values = np.loadtxt(modelopts["model_params_file"])
    p_lines = [ x.rstrip("\n") for x in open(modelopts["pairwise_params_file"],"r").readlines() ]

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
        pairwise_other_params.append(tuple(temp))

    pairs = np.array(pairs) 
    pairwise_param_assignment = np.array(pairwise_param_assignment)
    pairwise_type = np.array(pairwise_type)
    modelopts["pairs"] = pairs
    modelopts["pairwise_param_assignment"] = pairwise_param_assignment
    modelopts["pairwise_type"] = pairwise_type
    modelopts["pairwise_other_parameters"] = pairwise_other_params
    modelopts["model_param_values"] = model_param_values 
    modelopts["pairwise_params_file_location"] = modelopts["pairwise_params_file"]
    modelopts["model_params_file_location"] = modelopts["model_params_file"]
    modelopts["defaults"] = False

def _add_pair_opts(modelopts):
    """Add pairs or pairwise_params info to modelopts"""
    if (modelopts["pairs_file"] == None) and (modelopts["pairwise_params_file"] == None):
        raise IOError("Need to specify either pairs_file or pairwise_params_file")
    elif (modelopts["pairs_file"] != None):
        pairs = np.loadtxt("%s" % modelopts["pairs_file"],dtype=int)
        # For removing the unnecessary columns from the smog contact map output
        if len(pairs[0,:]) == 4:
            pairs = pairs[:,np.array([1,3])]
        modelopts["pairs"] = pairs
        modelopts["defaults"] = True
    else:
        _add_pairwise_params(modelopts)
        
#############################################################################
# Define Functions specific to a package (i.e. FRET, ddG_MC2004)
#############################################################################

def FRET_fitopts_load(item, value):
    """Check options specific to FRET-package"""
    bool_valid_check = ["true", "True", "1", "T", "t"]
    if item == "t_fit":
        value = int(value)
    elif item == "fret_pairs":
        value = [ int(x) for x in re.split(",\s+|\s+", value.strip("[ | ]"))]
        if (len(value) % 2) != 0:
            raise IOError("len(fret_pairs) should be even. Invalid input: %s " % value.__repr__())
        temp_value = np.zeros((len(value)/2,2)).astype(int)
        for i in range(len(value)):
            temp_value[i/2, i%2] = value[i]
        value = temp_value
    elif item == "spacing":
        value = float(value)
    elif item == "truncate_value":
        value = float(value)
    elif item == "y_shift":
        value = float(value)
    elif item == "fretdata":
        value = re.split(",\s+|\s+", value.strip("[ | ] | '")) 
    elif item == "prevent_zero":
        value = value in bool_valid_check
    return value

def FRET_fitopts_save(key, option, config):
    if key == "fret_pairs":
        collection = ""
        print option
        for i in range(np.shape(option)[0]):
            for j in range(np.shape(option)[1]):
                collection += "%d " % option[i,j]
        config.set("fitting",key,collection)
    if key == "fretdata":
        collection = ""
        for string in option:
            collection += "%s " % string
        config.set("fitting",key,collection)
    
def tmatrix_fitopts_load(item,value):
    """check option for transition matrix"""
    if item == "lag_step": ##number of frames per a lag time step
        value = int(value)
    else:
        value = FRET_fitopts_load(item, value)
    return value



