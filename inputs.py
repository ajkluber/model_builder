""" Check inputs for making a CoarseGrainedModel """

import numpy as np
import re
import os
import shutil
import ConfigParser

import models.CoarseGrainedModel as cg

#############################################################################
# Helper functions to load in models from .ini files
#############################################################################
def save_model(model,fitopts):
    name = model.name
    config = ConfigParser.SafeConfigParser(allow_no_value=True)
    config.add_section("model")
    config.add_section("fitting")
    modelkeys = ["name","bead_repr","disulfides","pairs_file",
                "pairwise_params_file_location","model_params_file_location",
                "defaults","cb_volume","n_native_pairs","contact_type",
                "backbone_param_vals","starting_gro","simple_disulfides",
                "verbose","using_sbm_gmx"]

    # Save fitting options that aren't None.
    for key in fitopts.iterkeys():
        if fitopts[key] not in [None,""]:
            if key == "include_dirs":
                temp = ""
                for dir in fitopts["include_dirs"]:
                    temp += "%s " % dir
                config.set("fitting",key,temp)
            else:
                config.set("fitting",key,str(fitopts[key]))
    
    # Save model options that aren't None.
    for key in modelkeys:
        value = getattr(model,key)
        if value not in [None,""]:
            if key in ["pairwise_params_file_location","model_params_file_location"]:
                config.set("model",key.split("_location")[0],str(value))
            else:
                config.set("model",key,str(value))

    if os.path.exists("%s.ini" % name):
        shutil.move("%s.ini" % name,"%s.1.ini" % name)

    with open("%s.ini" % name,"w") as cfgfile:
        config.write(cfgfile)

def load_model(name,dry_run=False):
    modelopts, fittingopts = load_config(name)
    modelopts["dry_run"] = dry_run
    model = cg.CoarseGrainedModel(**modelopts)
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

def new_models(names):
    """ Create new models with inputted options."""
    Models = []
    Fittingopts = []
    for name in names:
        modelopts, fittingopts = load_config(name)
        model = cg.CoarseGrainedModel(**modelopts)
        Models.append(model)
        Fittingopts.append(fittingopts)
    return Models,Fittingopts

def new_model_from_config(name):
    """ Create new models with inputted options."""
    modelopts, fittingopts = load_config(name)
    model = cg.CoarseGrainedModel(**modelopts)
    return model

def get_pairwise_params(pairwise_params_file,model_params_file):
    """ Grab pairwise_params from file. """
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
        pairwise_other_params.append(tuple(temp))

    pairs = np.array(pairs) 
    pairwise_param_assignment = np.array(pairwise_param_assignment)
    pairwise_type = np.array(pairwise_type)

    return pairs,pairwise_param_assignment,model_param_values,pairwise_type,pairwise_other_params

def load_config(name):
    """Parse options from <name>.ini file"""
    config = ConfigParser.SafeConfigParser(allow_no_value=True)
    config.read("%s.ini" % name)

    # Populate all fields as None 
    modelopts = _empty_model_opts()
    fittingopts = _empty_fitting_opts()

    print "Creating model according to %s.ini" % name
    print "Options not shown default to None"
    load_model_section(config.items("model"),modelopts)
    load_fitting_section(config,modelopts,fittingopts)
    _add_pair_opts(modelopts) 
    modelopts["pdb"] = "%s.pdb" % modelopts["name"]
    return modelopts,fittingopts

def load_model_section(modelitems,modelopts):
    """Parse [model] options from config .ini file"""
    print "Model options:"
    for item,value in modelitems:
        if value in [None,""]:
            pass
        else:
            print "  %-20s = %s" % (item,value)
            if item == "n_native_pairs":
                value = int(value)
            elif item == "epsilon_bar":
                value = float(value)
            elif item == "using_sbm_gmx":
                value = bool(value)
            elif item == "umbrella":
                value = bool(value)
            elif item == "simple_disulfides":
                value = bool(value)
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
            elif item == "name":
                if not os.path.exists("%s.pdb" % value):
                    raise IOError("%s.pdb file does not exist! Check config file inputs" % value)
            elif item == "starting_gro":
                if not os.path.exists("%s" % value):
                    raise IOError("%s file does not exist! Check config file inputs" % value)
            modelopts[item] = value

    if modelopts["bead_repr"] == "CACB" and (modelopts["cb_volume"] not in ["average","flavored"]):
        raise IOError("If bead_repr = CACB, then must set cb_volume to: average or flavored")

def load_fitting_section(config,modelopts,fittingopts):
    """Parse [fitting] options from config .ini file"""
    # special fitting checks is for package specific options
    # assigns based on keys, functions should be at end of file
    special_fitting_checks = {"FRET":FRET_fitopts_load, "tmatrix":tmatrix_fitopts_load}
        
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
                    value = bool(value)
                elif item == "constrain_avg_eps":
                    value = bool(value)
                elif item == "nonnative":
                    value = bool(value)
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
                    value = bool(value)
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
            "using_sbm_gmx","umbrella"] 
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
    if item == "t_fit":
        value = int(value)
    elif item == "fret_pairs":
        value = [ int(x) for x in re.split(",\s+|\s+", value.strip("[ | ]"))]
        if (len(value) % 2) != 0:
            raise IOError("len(fret_pairs) should be even. Invalid input: %s " % value.__repr__())
        temp_value = []
        holder = [0,0]
        for i in range(len(value)):
            if i%2 == 0:
                holder[0] = value[i]
            else:
                holder[1] = value[i]
                temp_value.append(holder)
        value = temp_value
    elif item == "spacing":
        value = float(value)
    elif item == "truncate_value":
        value = float(value)
    elif item == "y_shift":
        value = float(value)
    elif item == "fretdata":
        value = str(value)
    
    return value

def tmatrix_fitopts_load(item,value):
    """check option for transition matrix"""
    if item == "lag_step": ##number of frames per a lag time step
        value = int(value)
    else:
        value = FRET_fitopts_load(item, value)
    return value



