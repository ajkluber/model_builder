
def get_model_info(model):
    """ The string representation of all the model info."""
    model_info_string = "[ Path ]\n"
    model_info_string += "%s\n" % model.path
    model_info_string += "[ PDB ]\n"
    model_info_string += "%s\n" % model.pdb
    model_info_string += "[ Subdir ]\n"
    model_info_string += "%s\n" % model.subdir
    model_info_string += "[ Iteration ]\n"
    model_info_string += "%d\n" % model.iteration
    model_info_string += "[ Model_Code ]\n" 
    model_info_string += "%s\n" % model.model_code
    model_info_string += "[ Bead_Model ]\n" 
    model_info_string += "%s\n" % model.bead_repr
    model_info_string += "[ Backbone_params ]\n" 
    model_info_string += "  %5s       %5s       %5s\n" % ("Kb","Ka","Kd")
    model_info_string += "[ Backbone_param_vals ]\n" 
    model_info_string += "%10.2f%10.2f%10.2f\n" % \
        (model.backbone_param_vals["Kb"],model.backbone_param_vals["Ka"],model.backbone_param_vals["Kd"])
    model_info_string += "[ Disulfides ]\n"
    if model.disulfides == None:
        model_info_string += "%s\n" % None
    else:
        temp = ''
        for x in model.disulfides:
            temp += " %d " % x 
        model_info_string += "%s\n" % temp
    model_info_string += "[ Epsilon_Bar ]\n"
    model_info_string += "%s\n" % str(model.epsilon_bar)
    model_info_string += "[ Nonnative ]\n"
    model_info_string += "%s\n" % str(model.nonnative)
    model_info_string += "[ Pairwise_Params_File ]\n"
    model_info_string += "%s\n" % str(model.pairwise_params_file_location)
    model_info_string += "[ Model_Params_File ]\n"
    model_info_string += "%s\n" % str(model.model_params_file_location)
    model_info_string += "[ N_Native_Pairs ]\n"
    model_info_string += "%s\n" % str(model.n_native_pairs)
    model_info_string += "[ Contact_Type ]\n"
    model_info_string += "%s\n" % model.contact_type
    model_info_string += "[ Fitting_Data ]\n"
    model_info_string += "%s\n" % model.fitting_data
    model_info_string += "[ Fitting_Includes ]\n"
    for dir in model.fitting_includes:
        model_info_string += "%s" % str(dir)
    model_info_string += "\n"
    model_info_string += "[ Fitting_Solver ]\n"
    model_info_string += "%s\n" % model.fitting_solver
    model_info_string += "[ Fitting_AllowSwitch ]\n"
    model_info_string += "%s\n" % model.fitting_allowswitch
    model_info_string += "[ Fitting_Params_File ]\n"
    model_info_string += "%s\n" % model.fitting_params_file
    model_info_string += "[ Reference ]\n" 
    model_info_string += "None\n" 
    return model_info_string
