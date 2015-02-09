import os
import argparse
import ConfigParser
import glob

def convert_info_to_config(name):

    # Config file should at least have [model] section.
    # All options that aren't specified should be set to None and will be left
    # as empty when config file is written

    cfg = ConfigParser.SafeConfigParser(allow_no_value=True)
    cfg.add_section("model")
    cfg.add_section("fitting")
    cfg.set("model","name",name)

    info_file = open('%s/model.info' % name,'r')
    lasttime,action,task = open('%s/modelbuilder.log' % name,'r').readlines()[-1].split()
    last_comp_task = "%s %s" % (action,task)
    info = [ x.rstrip("\n") for x in info_file.readlines() ]
    options = { info[2*i].split()[1]:info[2*i + 1] for i in range(len(info)/2) }
    
    opts = ["Fitting_Data","Fitting_Includes","Fitting_Solver",
            "Iteration","Fitting_AllowSwitch","Fitting_Params_File",
            "Nonnative","Pairs_File","Pairwise_Params_File",
            "Model_Params_File","Epsilon_Bar",
            "Defaults","Disulfides",
            "N_Native_Pairs","Contact_Type","Model_Code"]
    for opt in opts:
        if not options.has_key(opt):
            options[opt] = "None"

    # [model] section
    # Required [model] fields are:
    #   1. pdb/name
    #   2. bead_repr
    #   3. contacts_file or pairwise_params_file
    #   
    # Optional fields:
    #   4. code 
    #   5. backbone_params 
    #   6. backbone_param_vals 
    #   7. disulfides
    #   8. contact_type
    #   9. n_native_contacts
    #   10. nonative 
    #   11. epsilon_bar
    cfg.set("model","bead_repr",options["Bead_Model"]) 

    if options["Disulfides"] != "None":
        cfg.set("model","disulfides",options["Disulfides"]) 

    if options["Pairwise_Params_File"] != "None":
        cfg.set("model","pairwise_params_file",options["Pairwise_Params_File"]) 

    if options["Pairwise_Params_File"] != "None":
        cfg.set("model","model_params_file",options["Model_Params_File"]) 
        cfg.set("model","defaults","False") 

    if options["Contact_Type"] != "None":
        cfg.set("model","contact_type",options["Contact_Type"]) 

    if options["N_Native_Pairs"] != "None":
        cfg.set("model","n_native_pairs",options["N_Native_Pairs"]) 

    if options["Epsilon_Bar"] != "None":
        cfg.set("model","epsilon_bar",options["Epsilon_Bar"]) 

    # [fitting] section
    # Fields for [fitting]
    #   1. data_type
    #   2. include_dirs
    #   3. solver
    #   4. allow_switch
    #   5. params_to_fit
    if options["Fitting_Data"] != "None":
        cfg.set("fitting","data_type",options["Fitting_Data"]) 

    if options["Fitting_Includes"] not in ["None",name]:
        cfg.set("fitting","include_dirs",options["Fitting_Includes"]) 

    if options["Fitting_Solver"] != "None":
        cfg.set("fitting","solver",options["Fitting_Solver"]) 

    if options["Iteration"] != "None":
        cfg.set("fitting","iteration",options["Iteration"]) 

    if options["Fitting_AllowSwitch"] != "None":
        cfg.set("fitting","allow_switch",options["Fitting_AllowSwitch"]) 

    if options["Fitting_Params_File"] != "None":
        cfg.set("fitting","parameters_to_fit",options["Fitting_Params_File"]) 

    if options["Nonnative"] != "None":
        cfg.set("fitting","nonnative",options["Nonnative"]) 

    cfg.set("fitting","last_completed_task",last_comp_task) 
    
    with open("%s.ini" % name,"w") as cfgfile:
        cfg.write(cfgfile)

if __name__ == "__main__":
    files = glob.glob("*/model.info")
    for file in files:
        name = file.split("/")[-2]
        convert_info_to_config(name)

