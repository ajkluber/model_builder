import numpy as np

import bonded_potentials as bond
import pdb_parser

############################################################################
# Helper function to get representation
############################################################################
def set_bonded_interactions(model):
    allowed = "(CA, CACB)"
    if model.bead_repr == "CA":
        set_CA_bonded_interactions(model)
    elif model.bead_repr== "CACB":
        set_CACB_bonded_interactions(model)
    else:
        raise IOError('model.bead_repr must be either %s.' % allowed)

    check_disulfides(model)
    get_atoms_string(model)
    get_bonds_string(model)
    get_angles_string(model)
    get_dihedrals_string(model)
    get_grofile(model)

############################################################################
# Calpha representation
############################################################################
def set_CA_bonded_interactions(model):
    """ Extract info from the Native.pdb for making index and top file """
    # Grab coordinates from the pdb file.
    model.cleanpdb = pdb_parser.get_clean_CA(model.pdb)
    model.cleanpdb_full = pdb_parser.get_clean_full(model.pdb)
    model.cleanpdb_full_noH = pdb_parser.get_clean_full_noH(model.pdb)
    pdb_info = pdb_parser.get_coords_atoms_residues(model.cleanpdb)

    model.atm_coords = pdb_info[0]
    model.atm_indxs = pdb_info[1]
    model.atm_types = pdb_info[2]
    model.res_indxs = pdb_info[3]
    model.res_types = pdb_info[4]
    indxs = model.atm_indxs
    coords = model.atm_coords

    model.n_residues = len(np.unique(np.array(model.res_indxs,copy=True)))
    model.n_atoms = len(model.atm_types)

    # Set bonded force field terms quantities.
    model.bond_indices = [[indxs[i],indxs[i+1]] for i in range(model.n_atoms-1)]
    model.angle_indices = [[indxs[i],indxs[i+1],indxs[i+2]] for i in range(model.n_atoms-2)]
    model.dihedral_indices = [[indxs[i],indxs[i+1],indxs[i+2],indxs[i+3]] for i in range(model.n_atoms-3)]

    model.bond_min = [ bond.distance(coords,i_idx-1,j_idx-1) for i_idx,j_idx in model.bond_indices ]
    model.angle_min = [ bond.angle(coords,i_idx-1,j_idx-1,k_idx-1) for i_idx,j_idx,k_idx in model.angle_indices ]
    model.dihedral_min = [ bond.dihedral(coords,i_idx-1,j_idx-1,k_idx-1,l_idx-1) for i_idx,j_idx,k_idx,l_idx in model.dihedral_indices ]

    if not hasattr(model,"bond_strengths"):
        model.bond_strengths = [ model.backbone_param_vals["Kb"] for i in range(len(model.bond_min)) ]
    if not hasattr(model,"angle_strengths"):
        model.angle_strengths = [ model.backbone_param_vals["Ka"] for i in range(len(model.angle_min)) ]
    if not hasattr(model,"dihedral_strengths"):
        model.dihedral_strengths = [ model.backbone_param_vals["Kd"] for i in range(len(model.dihedral_min)) ]

    # Using proper dihedrals for all C-alpha dihedrals.
    model.dihedral_type = [ 1 for i in range(len(model.dihedral_min)) ]

    # atomtypes category of topol.top. Sets default excluded volume of 0.4nm
    atomtypes_string = " [ atomtypes ]\n"
    atomtypes_string += " ;name  mass     charge   ptype c10       c12\n"
    atomtypes_string += " CA     1.000    0.000 A    0.000   %10.9e\n\n" % (0.4**12)
    model.atomtypes_string = atomtypes_string

    # Create index.ndx string for CA representation.
    ca_string = ''
    i = 1
    for ind in indxs: 
        if (i % 15) == 0:
            ca_string += '%4d \n' % ind
        else:
            ca_string += '%4d ' % ind
        i += 1
    ca_string += '\n'
    headings = ["System","Protein","Protein-H","C-alpha",\
                "Backbone","MainChain","MainChain+Cb","MainChain+H"]
    indexstring = ""
    for heading in headings:
        indexstring += "[ "+heading+" ]\n"
        indexstring += ca_string
    indexstring += '[ SideChain ]\n\n'
    indexstring += '[ SideChain-H ]\n\n'
    model.index_ndx = indexstring

############################################################################
# Calpha Cbeta representation
############################################################################
def set_CACB_bonded_interactions(model):
    """Extract info from the Native.pdb for making index and top file. 
    NOT DONE
    """
    # Grab coordinates from the pdb file.
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

    # atomtypes category of topol.top. Sets default excluded volume of ???
    atomtypes_string = " [ atomtypes ]\n"
    atomtypes_string += " ;name  mass     charge   ptype c10       c12\n"
    atomtypes_string += " CA     1.000    0.000 A    0.000   %10.9e\n" % (0.19**12)
    atomtypes_string += " CB     1.000    0.000 A    0.000   %10.9e\n\n" % (0.1**12) 
    model.atomtypes_string = atomtypes_string

    # Make index.ndx string
    create_CACB_index_ndx(model,CA_indxs,CB_indxs)

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
    """To keep torsional potential in same proportion to nonbonded potential"""
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
    Cheung 'Exploring the interplay between topology and secondary structural
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

############################################################################
# Functions to check that disulfides are reasonable
############################################################################
def check_disulfides(model):
    """ Check that specified disulfides are between cysteine and that 
        the corresonding pairs are within 0.8 nm. """
    coords = model.atm_coords
    residues = model.res_types
    if model.disulfides != None:
        if model.verbose:
            print "  Checking if disulfides are reasonable."
        for i in range(len(model.disulfides[::2])):
            i_idx = model.disulfides[2*i]
            j_idx = model.disulfides[2*i + 1]
            dist = bond.distance(coords,i_idx-1,j_idx-1)
            theta1 = bond.angle(coords,i_idx-2,i_idx-1,j_idx-1)
            theta2 = bond.angle(coords,i_idx-1,j_idx-1,j_idx-2)
            phi = bond.dihedral(coords,i_idx-2,i_idx-1,j_idx-1,j_idx-2)
            if (residues[i_idx-1] != "CYS") or (residues[j_idx-1] != "CYS"):
                print "WARNING! Specifying disulfide without cysteines: %s  %s" % \
                        (residues[i_idx-1]+str(i_idx), residues[j_idx-1]+str(j_idx))
            if dist > 0.8:
                print "WARNING! Specifying disulfide with separation greater than 0.8 nm."
            else:
                if model.verbose:
                    print "   %s %s separated by %.4f nm, Good." % \
                (residues[i_idx-1]+str(i_idx), residues[j_idx-1]+str(j_idx),dist)
                # Remove disulfide pair from model.pairs if it is there.
                new_pairs = []
                for pair in model.pairs:
                    if (pair[0] == i_idx) and (pair[1] == j_idx):
                        continue
                    else:
                        new_pairs.append(pair)
                model.pairs = np.array(new_pairs)
                model.n_pairs = len(model.pairs)

                model.exclusions.append([i_idx,j_idx])
                # Set cysteine bond distance, angles, and dihedral.
                model.bond_indices.append([i_idx,j_idx])
                model.bond_min.append(dist)
                model.bond_strengths.append(model.backbone_param_vals["Kb"])

                model.angle_indices.append([i_idx-1,i_idx,j_idx])
                model.angle_indices.append([i_idx,j_idx,j_idx-1])
                model.angle_min.append(theta1)
                model.angle_min.append(theta2)
                model.angle_strengths.append(model.backbone_param_vals["Ka"])
                model.angle_strengths.append(model.backbone_param_vals["Ka"])

                model.dihedral_indices.append([i_idx-1,i_idx,j_idx,j_idx-1])
                model.dihedral_min.append(phi)
                model.dihedral_strengths.append(model.backbone_param_vals["Kd"])
    else:
        if model.verbose:
            print "  No disulfides to check."

############################################################################
# Functions to create the bonded sections of the topol.top and .gro file
############################################################################
def get_atoms_string(model):
    """ Generate the [ atoms ] string."""
    atoms_string = " [ atoms ]\n"
    atoms_string += " ;nr  type  resnr residue atom  cgnr charge  mass\n"
    for j in range(len(model.atm_indxs)):
        atmnum = model.atm_indxs[j]
        atmtype = model.atm_types[j]
        resnum = model.res_indxs[j]
        restype = model.res_types[j]
        atoms_string += " %5d%4s%8d%5s%4s%8d%8.3f%8.3f\n" % \
                    (atmnum,atmtype,resnum,restype,atmtype,atmnum,0.0,1.0)
    model.atoms_string = atoms_string

def get_bonds_string(model):
    """ Generate the [ bonds ] string."""
    bonds_string = " [ bonds ]\n"
    bonds_string += " ; ai aj func r0(nm) Kb\n"
    for k in range(len(model.bond_min)):
        i_idx = model.bond_indices[k][0]
        j_idx = model.bond_indices[k][1]
        dist = model.bond_min[k]
        kb = model.bond_strengths[k]
        bonds_string += "%6d %6d%2d%18.9e%18.9e\n" %  \
                      (i_idx,j_idx,1,dist,kb)
    model.bonds_string = bonds_string

def get_angles_string(model):
    """ Generate the [ angles ] string."""
    angles_string = " [ angles ]\n"
    angles_string += " ; ai  aj  ak  func  th0(deg)   Ka\n"
    for n in range(len(model.angle_min)):
        i_idx = model.angle_indices[n][0]
        j_idx = model.angle_indices[n][1]
        k_idx = model.angle_indices[n][2]
        theta = model.angle_min[n]
        ka = model.angle_strengths[n]
        angles_string += "%6d %6d %6d%2d%18.9e%18.9e\n" %  \
                      (i_idx,j_idx,k_idx,1,theta,ka)
    model.angles_string = angles_string

def get_dihedrals_string(model):
    """ Generate the [ dihedrals ] string."""
    dihedrals_string = " [ dihedrals ]\n"
    dihedrals_string += " ; ai  aj  ak al  func  phi0(deg)   Kd mult\n"
    dihedrals_ndx = '[ dihedrals ]\n'
    for n in range(len(model.dihedral_min)):
        i_idx = model.dihedral_indices[n][0]
        j_idx = model.dihedral_indices[n][1]
        k_idx = model.dihedral_indices[n][2]
        l_idx = model.dihedral_indices[n][3]
        dih_type = model.dihedral_type[n]
        phi = model.dihedral_min[n]
        kd = model.dihedral_strengths[n]
        if dih_type == 1:
            dihedrals_string += "%6d %6d %6d %6d%2d%18.9e%18.9e%2d\n" %  \
                          (i_idx,j_idx,k_idx,l_idx,dih_type,phi,kd,1)
            dihedrals_string += "%6d %6d %6d %6d%2d%18.9e%18.9e%2d\n" %  \
                          (i_idx,j_idx,k_idx,l_idx,dih_type,3.*phi,kd/2.,3)
        elif dih_type == 2:
            dihedrals_string += "%6d %6d %6d %6d%2d%18.9e%18.9e\n" %  \
                          (i_idx,j_idx,k_idx,l_idx,dih_type,phi,kd)
        dihedrals_ndx += '%4d %4d %4d %4d\n' % \
                            (i_idx,j_idx,k_idx,l_idx)
    model.dihedrals_string = dihedrals_string
    model.dihedrals_ndx = dihedrals_ndx

def get_grofile(model):
    """ Get the .gro string """
    gro_string = " Structure-based Gro file\n"
    gro_string += "%12d\n" % len(model.atm_types)
    for i in range(len(model.atm_indxs)):
        gro_string += "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" % \
            (model.res_indxs[i],model.res_types[i],
            model.atm_types[i],model.atm_indxs[i],
            model.atm_coords[i][0],model.atm_coords[i][1],model.atm_coords[i][2])
    gro_string += "   %-25.16f%-25.16f%-25.16f" % (50.0,50.0,50.0)

    model.grofile = gro_string

