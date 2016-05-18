import warnings

def interaction_exists_warning(pot):
    warnings.warn("Interaction already exists! skipping: {}".format(pot.describe()))

def default_sbm_parameters_warning():
    warnings.warn("Using default SBM parameters")

def default_sbm_potentials_warning():
    warnings.warn("Using default SBM parameters")

def missing_reference_warning():
    warnings.warn("Need to set reference structure model.set_reference()")
