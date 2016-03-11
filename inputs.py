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
    modelopts, fittingopts = load_config(name)
    modelopts["dry_run"] = dry_run
    topology = mdtraj.load(modelopts["topology"]).topology
    model = SBM(topology, bead_repr=modelopts["bead_repr"])
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
    check_exists = ["topology", "starting_gro", "name"]
    
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
            "using_sbm_gmx","umbrella", "topology"] 
    modelopts = { opt:None for opt in opts }
    return modelopts


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



