""" Check inputs for making a CoarseGrainedModel """

import numpy as np
import re
import os
import shutil
from glob import glob
import ConfigParser

import models.smog_AA_model as sg

#############################################################################
# Helper functions to load in models from .ini files
#############################################################################
def save_model(model,fitopts):
    name = model.name
    config = ConfigParser.SafeConfigParser(allow_no_value=True)
    config.add_section("model")
    config.add_section("fitting")
    modelkeys = ["name","bead_repr","long_pairs_file",
                "short_pairs_file", "long_pairs_file_location",
                "contact_type",
                "starting_gro",
                "verbose"]

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
            if key in ["long_pairs_file_location"]:
                config.set("model",key.split("_location")[0],str(value))
            else:
                config.set("model",key,str(value))

    if os.path.exists("%s.ini" % name):
        shutil.move("%s.ini" % name,"%s.1.ini" % name)

    with open("%s.ini" % name,"w") as cfgfile:
        config.write(cfgfile)

    model_path = model.path

#    smog_files = glob('{0}/{1}/smog_files/*'.format(model_path,name))
#    for item in smog_files:                  
#        shutil.copy(item,'.')

def load_model(name,dry_run=False):
    modelopts, fittingopts = load_config(name)
    modelopts["dry_run"] = dry_run
    model = sg.smog_AA_model(**modelopts)
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
        model = sg.smog_AA_model(**modelopts)
        Models.append(model)
        Fittingopts.append(fittingopts)
    return Models,Fittingopts

def new_model_from_config(name):
    """ Create new models with inputted options."""
    modelopts, fittingopts = load_config(name)
    model = sg.smog_AA_model(**modelopts)
    return model

def load_config(name):
    """Parse options from <name>.ini file"""
    if not os.path.exists("%s.ini" % name):
        raise IOError("%s.ini doesn't exist!" % name)
    config = ConfigParser.SafeConfigParser(allow_no_value=True)
    config.read("%s.ini" % name)

    # Populate all fields as None 
    modelopts = _empty_model_opts()
    fittingopts = _empty_fitting_opts()

    print "Creating model according to %s.ini" % name
    print "Options not shown default to None"
    load_model_section(config.items("model"),modelopts)
    load_fitting_section(config,modelopts,fittingopts)
    modelopts["pdb"] = "%s.pdb" % modelopts["name"]
    return modelopts,fittingopts

def load_model_section(modelitems,modelopts):
    """Parse [model] options from config .ini file"""
    print "Model options:"
    bool_valid_check = ["true", "True", "1", "T", "t"] #for checking boolean values
    for item,value in modelitems:
        if value in [None,""]:
            pass
        else:
            print "  %-20s = %s" % (item,value)
            if item == "epsilon_bar":
                value = float(value)
            elif item == "umbrella":
                value = value in bool_valid_check
            elif item.endswith("_file"):
                if not os.path.exists(value):
                    raise IOError("%s file does not exist! Check config file inputs" % value)
            elif item == "name":
                if not os.path.exists("%s.pdb" % value):
                    raise IOError("%s.pdb file does not exist! Check config file inputs" % value)
            elif item == "starting_gro":
                if not os.path.exists("%s" % value):
                    raise IOError("%s file does not exist! Check config file inputs" % value)
            elif item == "initial_t_array":
                value = list(value.split())
                value = [int(i) for i in value]
        
            modelopts[item] = value
    
def load_fitting_section(config,modelopts,fittingopts):
    """Parse [fitting] options from config .ini file"""
    # special fitting checks is for package specific options
    # assigns based on keys, functions should be at end of file
    bool_valid_check = ["true", "True", "1", "T", "t"]
#    queue_check = ["serial","commons","parallel","bigmem"]
    
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
                    value = value in bool_valid_check
                elif item == "constrain_avg_eps":
                    value = value in bool_valid_check
                elif item == "nonnative":
                    value = value in bool_valid_check
                elif item in ["equil_walltime","walltime"]:
                    if re.match("\d\d:\d\d:\d\d",value) == None:
                        raise IOError(" %s must be a time in format HH:MM:SS" % item)
                elif item == "queue":
                    value = value 
                elif item == "parameters_to_fit":
                    if not os.path.exists(value):
                        raise IOError("%s file does not exist! Check config file inputs" % value)
                    else:
                        modelopts["fitting_params"] = np.loadtxt(value,dtype=int)
                elif item == "cutoffs":
                    value = [ float(x) for x in re.split(",\s+|\s+", value.strip("[ | ]"))] 
                elif item == "simplify_lambdas":
                    value = value in bool_valid_check

                fittingopts[item] = value
                        
#############################################################################
# Internal functions to load in models from .ini files
#############################################################################
def _empty_fitting_opts():
    """Fitting options to check for"""
    opts = ["data_type","include_dirs","solver",
            "iteration","n_processors","queue"
            "allow_switch","parameters_to_fit",
            "nonnative","last_completed_task","cutoffs","simplify_lambdas"]         
    fittingopts = { opt:None for opt in opts }
    fittingopts["simplify_lambdas"] = False
    return fittingopts

def _empty_model_opts():
    """Model options to check for"""
    opts = ["name","bead_repr","long_pairs_file",
            "short_pairs_file","long_pairs_file_location","n_residues",
            "contact_type","starting_gro","epsilon_bar",
            "epsilon_bar","n_native_pairs",
            "n_short_native_pairs","pairs","long_pairs",
            "short_pairs","pairwise_other_parameters",
            "pairwise_param_assignment","n_long_native_pairs",
            "n_processors","n_pairs",
            "dihedrals","dihedrals_params","angles",
            "angles_params","bonds","bonds_params",
            "pairwise_type","verbose","dry_run",
            "umbrella", "pair_eps","pair_V",
            "long_pair_eps","long_pair_V",
            "n_fitting_params","model_param_values",
            "long_model_param_values","qref",
            "initial_t_array"
            ] 
    modelopts = { opt:None for opt in opts }
    return modelopts




