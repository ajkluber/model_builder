import numpy as np

import bonded_potentials as bond
import pdb_parser
import residue_properties as rp

############################################################################
# Calpha Cbeta representation
############################################################################
def set_CACB_bonded_interactions(model):
    """Extract info from the Native.pdb for making index and top file. 
    """
    # Grab useful info from the pdb file.
    model.cleanpdb = pdb_parser.get_clean_CA_center_of_mass_CB(model.pdb)
    model.cleanpdb_full = pdb_parser.get_clean_full(model.pdb)
    model.cleanpdb_full_noH = pdb_parser.get_clean_full_noH(model.pdb)
    pdb_info = pdb_parser.get_coords_atoms_residues(model.cleanpdb)

    model.atm_coords = pdb_info[0]
    model.atm_indxs = pdb_info[1]
    model.atm_types = pdb_info[2]
    model.res_indxs = pdb_info[3]
    model.res_types = pdb_info[4]
    model.res_indxs_unique = np.unique(np.array(model.res_indxs,copy=True))
    model.res_types_unique = pdb_info[5]
    indxs = model.atm_indxs
    coords = model.atm_coords


    model.n_residues = len(model.res_indxs_unique)
    model.n_atoms = len(model.atm_types)
    CA_indxs = model.atm_indxs[model.atm_types == "CA"]
    CB_indxs = model.atm_indxs[model.atm_types == "CB"]
    model.CA_indxs = CA_indxs
    model.CB_indxs = CB_indxs

    # Collect the indices for atoms in bond, angle, and dihedral interactions.
    create_CACB_bonds(model,CA_indxs,CB_indxs,coords)
    create_CACB_angles(model,CA_indxs,CB_indxs,coords)
    create_CACB_dihedrals_improper(model,CA_indxs,CB_indxs,coords)
    create_CACB_exclusions(model)

    # atomtypes category of topol.top. Sets default excluded volume to 0.2 nm
    ca_size = 0.28
    #model.cb_volume = "flavored" # flavored or average
    residue_custom = {}

    if model.cb_volume.endswith(".dat"):
        residue_names, residue_radii = np.loadtxt("%s.dat" % model.cb_volume, unpack=True)
        for i in range(len(residue_radii)):
            residue_custom[residue_names[i]] = residue_radii[i] # prepares custom radii dictionary
        
    model.res_types_abbrev = []
    model.atm_names = []
    model.atm_radii = []
    
    for i in range(len(model.atm_indxs)):
        model.res_types_abbrev.append(rp.residue_code[model.res_types[i]])
        model.atm_names.append(model.atm_types[i] + model.res_types_abbrev[i])
        if  model.atm_types[i] == "CB":
            if model.cb_volume == "flavored":
                model.atm_radii.append(rp.residue_cacb_effective_interaction[model.res_types[i]])
            elif model.cb_volume == "average":
                model.atm_radii.append(rp.residue_cacb_effective_interaction["AVERAGE"])
            elif model.cb_volume.endswith(".dat"):
                model.atm_radii.append(float(residue_custom[model.res_types[i]]))
            else:
                raise IOError('cb_volume must be flavored, average, or custom .dat input file')
        else:
            model.atm_radii.append(ca_size) # for all CA* variant
    
    model.atm_names_dups , model.atm_names_dups_ind = np.unique(model.atm_names,return_index=True)

    atomtypes_string = " [ atomtypes ]\n"
    atomtypes_string += " ;name  mass     charge   ptype c10       c12\n"
    for i in model.atm_names_dups_ind[:-1]:
        atomtypes_string += (model.atm_names[i] + "     1.000    0.000 A    0.000   %10.9e\n" % (model.atm_radii[i]**12))
    atomtypes_string += (model.atm_names[model.atm_names_dups_ind[-1]] + "     1.000    0.000 A    0.000   %10.9e\n\n" % (model.atm_radii[model.atm_names_dups_ind[-1]]**12))
    model.atomtypes_string = atomtypes_string
    # Make index.ndx string
    create_CACB_index_ndx(model,CA_indxs,CB_indxs)
    return model.atm_names

def create_CACB_bonds(model,CA_indxs,CB_indxs,coords):
    # Set bonded force field terms quantities.
    model.bond_min = [] 
    model.bond_indices = []
    for i in range(model.n_residues-1):
        # First bond all the c-alphas together.
        model.bond_indices.append([CA_indxs[i],CA_indxs[i+1]])
        bond_dist = bond.distance(coords,CA_indxs[i]-1,CA_indxs[i+1]-1)
        model.bond_min.append(bond_dist)
    sub = 0
    for i in range(model.n_residues):
        # Then bond all c-alphas to their c-beta.
        if model.res_types_unique[i] == "GLY":
            # Skip glycine
            sub += 1
        else:
            model.bond_indices.append([CA_indxs[i],CB_indxs[i-sub]])
            bond_dist = bond.distance(coords,CA_indxs[i]-1,CB_indxs[i-sub]-1)
            model.bond_min.append(bond_dist)
    if not hasattr(model,"bond_strengths"):
        model.bond_strengths = [ model.backbone_param_vals["Kb"] for i in range(len(model.bond_min)) ]

def create_CACB_angles(model,CA_indxs,CB_indxs,coords):
    # Set angle terms
    model.angle_indices = []
    model.angle_min = []
    for i in range(model.n_residues-2):
        # First set angles for c-alphas.
        model.angle_indices.append([CA_indxs[i],CA_indxs[i+1],CA_indxs[i+2]])
        angle = bond.angle(coords,CA_indxs[i]-1,CA_indxs[i+1]-1,CA_indxs[i+2]-1)
        model.angle_min.append(angle)
    sub = 0
    for i in range(model.n_residues):
        # Then set angles between c-alphas and c-betas. 
        if model.res_types_unique[i] == "GLY":
            # Count how many glycines 
            sub += 1 
        else:
            if i == 0:
                # Account for N-terminus
                model.angle_indices.append([CB_indxs[i],CA_indxs[i],CA_indxs[i+1]])
                angle = bond.angle(coords,CB_indxs[i]-1,CA_indxs[i]-1,CA_indxs[i+1]-1)
                model.angle_min.append(angle)
            elif i == (model.n_residues - 1):
                # Account for C-terminus
                model.angle_indices.append([CA_indxs[i-1],CA_indxs[i],CB_indxs[i-sub]])
                angle = bond.angle(coords,CA_indxs[i-1]-1,CA_indxs[i]-1,CB_indxs[i-sub]-1)
                model.angle_min.append(angle)
            else:
                # Account for the middle.
                model.angle_indices.append([CA_indxs[i-1],CA_indxs[i],CB_indxs[i-sub]])
                model.angle_indices.append([CB_indxs[i-sub],CA_indxs[i],CA_indxs[i+1]])
                angle1 = bond.angle(coords,CA_indxs[i-1]-1,CA_indxs[i]-1,CB_indxs[i-sub]-1)
                angle2 = bond.angle(coords,CB_indxs[i-sub]-1,CA_indxs[i]-1,CA_indxs[i+1]-1)
                model.angle_min.append(angle1)
                model.angle_min.append(angle2)

    if not hasattr(model,"angle_strengths"):
        model.angle_strengths = [ model.backbone_param_vals["Ka"] for i in range(len(model.angle_min)) ]

def create_CACB_dihedrals_improper(model,CA_indxs,CB_indxs,coords):
    """Create chirality enforcing dihedrals using improper dihedrals."""
    # Set dihedral terms
    model.dihedral_indices = []
    model.dihedral_min = []
    model.dihedral_type = []
    model.dihedral_strengths = [] 
    for i in range(model.n_residues-3):
        # First set dihedrals for c-alphas.
        # Using proper dihedrals for all C-alpha dihedrals.
        model.dihedral_indices.append([CA_indxs[i],CA_indxs[i+1],CA_indxs[i+2],CA_indxs[i+3]])
        dihedral = bond.dihedral(coords,CA_indxs[i]-1,CA_indxs[i+1]-1,CA_indxs[i+2]-1,CA_indxs[i+3]-1)
        model.dihedral_min.append(dihedral)
        model.dihedral_type.append(1)
        model.dihedral_strengths.append(model.backbone_param_vals["Kd"])
    sub = 0
    for i in range(model.n_residues):
        # Then set dihedrals between c-alphas and c-betas. 
        if model.res_types_unique[i] == "GLY":
            # Counted skipped glycines to keep C-beta indices correct.
            sub += 1
        else:
            # Using improper (harmonic) dihedrals to enforce chirality on C-betas.
            if (i > 0) and (i < model.n_residues - 1):
                model.dihedral_indices.append([CA_indxs[i],CA_indxs[i-1],CA_indxs[i+1],CB_indxs[i-sub]])
                dih = bond.improper_dihedral(coords,CA_indxs[i]-1,CA_indxs[i-1]-1,CA_indxs[i+1]-1,CB_indxs[i-sub]-1)
                model.dihedral_min.append(dih)
                model.dihedral_type.append(2)    
                model.dihedral_strengths.append(model.backbone_param_vals["Ka"])

def create_CACB_dihedrals_proper(model,CA_indxs,CB_indxs,coords):
    """Create chirality enforcing dihedrals using proper dihedrals. NEVER USED."""

    # Set dihedral terms
    model.dihedral_indices = []
    model.dihedral_min = []
    for i in range(model.n_residues-3):
        # First set dihedrals for c-alphas.
        model.dihedral_indices.append([CA_indxs[i],CA_indxs[i+1],CA_indxs[i+2],CA_indxs[i+3]])
        dihedral = bond.dihedral(coords,CA_indxs[i]-1,CA_indxs[i+1]-1,CA_indxs[i+2]-1,CA_indxs[i+3]-1)
        model.dihedral_min.append(dihedral)
    sub = 0
    for i in range(model.n_residues-1):
        # Then set dihedrals between c-alphas and c-betas. 
        if model.res_types_unique[i] == "GLY":
            # Counted skipped glycines to keep C-beta indices correct.
            sub += 1
        else:
            if (i >= 2) and (i <= (model.n_residues - 3)):
                # Account for the middle
                model.dihedral_indices.append([CB_indxs[i-sub],CA_indxs[i],CA_indxs[i+1],CA_indxs[i+2]])
                model.dihedral_indices.append([CA_indxs[i-2],CA_indxs[i-1],CA_indxs[i],CB_indxs[i-sub]])
                dih_A = bond.dihedral(coords,CB_indxs[i-sub]-1,CA_indxs[i]-1,CA_indxs[i+1]-1,CA_indxs[i+2]-1)
                dih_B = bond.dihedral(coords,CA_indxs[i-2]-1,CA_indxs[i-1]-1,CA_indxs[i]-1,CB_indxs[i-sub]-1)
                model.dihedral_min.append(dih_A)
                model.dihedral_min.append(dih_B)
                # Check if glycine is next residue.
                if model.res_types_unique[i+1] != "GLY":
                    model.dihedral_indices.append([CB_indxs[i-sub],CA_indxs[i],CA_indxs[i+1],CB_indxs[i+1-sub]])
                    dih_C = bond.dihedral(coords,CB_indxs[i-sub]-1,CA_indxs[i]-1,CA_indxs[i+1]-1,CB_indxs[i+1-sub]-1)
                    model.dihedral_min.append(dih_C)
            else:
                # Note: Check if glycine is next residue.
                if i < 2:
                    # Account for N-terminus
                    model.dihedral_indices.append([CB_indxs[i-sub],CA_indxs[i],CA_indxs[i+1],CA_indxs[i+2]])
                    dih_A = bond.dihedral(coords,CB_indxs[i-sub]-1,CA_indxs[i]-1,CA_indxs[i+1]-1,CA_indxs[i+2]-1)
                    model.dihedral_min.append(dih_A)
                    # Check if glycine is next residue.
                    if model.res_types_unique[i+1] != "GLY":
                        model.dihedral_indices.append([CB_indxs[i-sub],CA_indxs[i],CA_indxs[i+1],CB_indxs[i+1-sub]])
                        dih_C = bond.dihedral(coords,CB_indxs[i-sub]-1,CA_indxs[i]-1,CA_indxs[i+1]-1,CB_indxs[i+1-sub]-1)
                        model.dihedral_min.append(dih_C)
                elif i == (model.n_residues - 2):
                    # Account for C-terminus
                    if model.res_types_unique[i+1] != "GLY":
                        model.dihedral_indices.append([CB_indxs[i-sub],CA_indxs[i],CA_indxs[i+1],CB_indxs[i+1-sub]])
                        dihedral2 = bond.dihedral(coords,CB_indxs[i-sub]-1,CA_indxs[i]-1,CA_indxs[i+1]-1,CB_indxs[i+1-sub]-1)
                        model.dihedral_min.append(dihedral2)
    # Using proper dihedrals for all C-alpha dihedrals.
    model.dihedral_type = [ 1 for i in range(len(model.dihedral_min)) ]

def scale_dihedral_strengths(model,CA_indxs):
    """To keep torsional potential in same proportion to nonbonded potential. NEVER USED"""
    model.dihedral_strengths = np.zeros(len(model.dihedral_indices),float)
    for i in range(model.n_residues-1):
        # Count the dihedrals that subsequent C-alphas are a part of.
        ca1 = CA_indxs[i]
        ca2 = CA_indxs[i+1]
        dih_count = 0
        dih_indxs = []
        for n in range(len(model.dihedral_indices)):
            if (ca1 == model.dihedral_indices[n][1]) and (ca2 == model.dihedral_indices[n][2]):
                dih_indxs.append(n)
                dih_count += 1
        # Reduce the dihedral strength for those dihedrals by their 
        # number
        model.dihedral_strengths[dih_indxs] = model.backbone_param_vals["Kd"]/float(dih_count)

def create_CACB_exclusions(model):
    """Create list of exclusions for the CACB model

    Rules for exclusions are based off the following reference:
    Cheung,et.al. 'Exploring the interplay between topology and secondary structural
        formation in the protein folding problem'. J.Chem.Phys. B. 2003.
    """

    if not hasattr(model,"exclusions"):
        model.exclusions = []

    # Exclude neighbors closer in sequence than:
    #   |i - j| < 4 for CA_i CA_i pairs
    #   |i - j| < 2 for CB_i CB_j pairs
    #   |i - j| < 2 for CA_i CB_j pairs
    cutAA = 4
    cutAB = 2
    cutBB = 2

    # Loop over all possible pairs
    for i in range(model.n_atoms):
        resid1 = model.res_indxs[i]
        indx1 = model.atm_indxs[i]
        type1 = model.atm_indxs[i]
        for j in range(model.n_atoms):
            resid2 = model.res_indxs[j]
            indx2 = model.atm_indxs[j]
            type2 = model.atm_indxs[j]
            excl_pair = [indx1,indx2]
            excl_pairr = [indx2,indx1]
            if i == j:
                continue

            # Add to exclusions if pair is closer in sequence
            # then allowed by the corresponding cutoff.
            if (type1 == "CA") and (type2 == "CA"):
                if abs(resid1 - resid2) < cutAA:
                    if (excl_pair not in model.exclusions) and \
                       (excl_pairr not in model.exclusions):
                        model.exclusions.append(excl_pair)
            elif (type1 == "CB") and (type2 == "CB"):
                if abs(resid1 - resid2) < cutBB:
                    if (excl_pair not in model.exclusions) and \
                       (excl_pairr not in model.exclusions):
                        model.exclusions.append(excl_pair)
            else:
                if abs(resid1 - resid2) < cutBB:
                    if (excl_pair not in model.exclusions) and \
                       (excl_pairr not in model.exclusions):
                        model.exclusions.append(excl_pair)


def create_CACB_index_ndx(model,CA_indxs,CB_indxs):
    """Record the indices for C-alpha, C-beta, and all atoms"""
    all_string = ''
    i = 1
    for ind in model.atm_indxs: 
        if (i % 15) == 0:
            all_string += '%4d \n' % ind
        else:
            all_string += '%4d ' % ind
        i += 1
    all_string += '\n'

    ca_string = ''
    i = 1
    for ind in CA_indxs: 
        if (i % 15) == 0:
            ca_string += '%4d \n' % ind
        else:
            ca_string += '%4d ' % ind
        i += 1
    ca_string += '\n'
    cb_string = ''
    i = 1
    for ind in CB_indxs: 
        if (i % 15) == 0:
            cb_string += '%4d \n' % ind
        else:
            cb_string += '%4d ' % ind
        i += 1
    cb_string += '\n'

    indexstring = "[ System ]\n"
    indexstring += all_string
    indexstring += "[ Protein ]\n"
    indexstring += all_string
    indexstring += "[ C-alpha ]\n"
    indexstring += ca_string
    indexstring += "[ Backbone ]\n"
    indexstring += ca_string
    indexstring += "[ C-beta ]\n"
    indexstring += cb_string
    indexstring += "[ Sidechain ]\n"
    indexstring += cb_string
    
    model.index_ndx = indexstring

