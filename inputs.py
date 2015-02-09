""" Check inputs for making a CoarseGrainedModel """

import numpy as np
import os
import shutil
import ConfigParser

import models.CoarseGrainedModel as cg
import convert_info_to_config as cvt

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
                "defaults","epsilon_bar",
                "n_native_pairs","contact_type",
                "verbose"]

    # Populate fields
    for key in fitopts.iterkeys():
        if fitopts[key] not in [None,""]:
            if key == "include_dirs":
                temp = ""
                for dir in fitopts["include_dirs"]:
                    temp += "%s " % dir
                config.set("fitting",key,temp)
            else:
                config.set("fitting",key,str(fitopts[key]))
    
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
    if not os.path.exists("%s.ini" % name):
        cvt.convert_info_to_config(name)
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
        pairwise_other_params.append(temp)

    pairs = np.array(pairs) 
    pairwise_param_assignment = np.array(pairwise_param_assignment)
    pairwise_type = np.array(pairwise_type)

    return pairs,pairwise_param_assignment,model_param_values,pairwise_type,pairwise_other_params

def load_config(name):
    """Load options from <name>.ini file"""
    config = ConfigParser.SafeConfigParser(allow_no_value=True)
    config.read("%s.ini" % name)

    # Populate all fields as None 
    modelopts = _empty_model_opts()
    fittingopts = _empty_fitting_opts()

    modelopts["pdb"] = "%s.pdb" % name
    print "Creating model according to %s.ini" % name
    print "Options not shown default to None"
    load_model_section(config,modelopts)
    load_fitting_section(config,modelopts,fittingopts)
    _add_pair_opts(modelopts) 
    return modelopts,fittingopts

def load_model_section(config,modelopts):
    print "Model options:"
    for item,value in config.items("model"):
        if value in [None,""]:
            pass
        else:
            print "  %-20s = %s" % (item,value)
            if item == "n_native_pairs":
                value = int(value)
            elif item == "epsilon_bar":
                value = float(value)
            elif item == "disulfides":
                import re
                value = [ int(x) for x in re.split(",\s+|\s+", value.strip("[ | ]"))]
                if (len(value) % 2) != 0:
                    raise IOError("len(disulfides) should be even. Invalid input: %s " % value.__repr__())
            elif item.endswith("_file"):
                if not os.path.exists(value):
                    raise IOError("%s file does not exist! Check config file inputs" % value)
            modelopts[item] = value

def load_fitting_section(config,modelopts,fittingopts):
    if config.has_section("fitting"):
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
                elif item == "nonnative":
                    value = bool(value)
                elif item == "parameters_to_fit":
                    if not os.path.exists(value):
                        raise IOError("%s file does not exist! Check config file inputs" % value)
                    else:
                        modelopts["fitting_params"] = np.loadtxt(value,dtype=int)
                ##specific to FRET-package
                elif item == "t_fit":
                    value = int(value)
                elif item == "fret_pairs":
                    import re
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
                    fittingopts[item] = temp_value
                elif item == "spacing":
                    value = float(value)
                    
#############################################################################
# Internal functions to load in models from .ini files
#############################################################################
def _empty_fitting_opts():
    opts = ["data_type","include_dirs","solver",
            "iteration","allow_switch","parameters_to_fit",
            "nonnative","last_completed_task"]         
    fittingopts = { opt:None for opt in opts }
    return fittingopts

def _empty_model_opts():
    opts = ["pairs_file","pairwise_params_file",
            "model_params_file","epsilon_bar",
            "defaults","bead_repr","disulfides",
            "n_native_pairs","contact_type","model_code",
            "pairs","pairwise_other_parameters",
            "pairwise_param_assignment",
            "pairwise_type","verbose","dry_run"]         
    modelopts = { opt:None for opt in opts }
    return modelopts


def _add_pairwise_params(modelopts):
    """Grab pairwise_params from file"""

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
        pairwise_other_params.append(temp)

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
    if (modelopts["pairs_file"] == None) and (modelopts["pairwise_params_file"] == None):
        raise IOError("Need to specify either pairs_file or pairwise_params_file")
    elif (modelopts["pairs_file"] != None):
        pairs = np.loadtxt("%s" % modelopts["pairs_file"],dtype=int)
        ##For removing the unnecessary columns from the smog contact map output
        if len(pairs[0,:]) == 4:
            pairs = pairs[:,np.array([1,3])]
        modelopts["pairs"] = pairs
        modelopts["defaults"] = True
    else:
        _add_pairwise_params(modelopts)

