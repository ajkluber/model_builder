import numpy as np

import bonded_potentials as bond
import pdb_parser
import cacb_repr

############################################################################
# Helper function to get representation
############################################################################
def set_bonded_interactions(model):
    """Create bonded section of the topol.top file.
    """
    allowed = "(CA, CACB)"
    if model.bead_repr == "CA":
        set_CA_bonded_interactions(model)
    elif model.bead_repr== "CACB":
        cacb_repr.set_CACB_bonded_interactions(model)
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
    model.atm_names = model.atm_types
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

                if not model.simple_disulfides:
                    print "not simple disulfides"
                    model.angle_indices.append([i_idx-1,i_idx,j_idx])
                    model.angle_indices.append([i_idx,j_idx,j_idx-1])
                    model.angle_min.append(theta1)
                    model.angle_min.append(theta2)
                    model.angle_strengths.append(model.backbone_param_vals["Ka"])
                    model.angle_strengths.append(model.backbone_param_vals["Ka"])

                    model.dihedral_indices.append([i_idx-1,i_idx,j_idx,j_idx-1])
                    model.dihedral_min.append(phi)
                    model.dihedral_type.append(1)
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
        atmtype = model.atm_names[j] # changed atm_types to atm_names
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
    gro_string += "%12d\n" % len(model.atm_names) # changed atm_types to atm_names
    for i in range(len(model.atm_indxs)):
        gro_string += "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" % \
            (model.res_indxs[i],model.res_types[i],
            model.atm_names[i],model.atm_indxs[i], # changed atm_types to atm_names
            model.atm_coords[i][0],model.atm_coords[i][1],model.atm_coords[i][2])
    gro_string += "   %-25.16f%-25.16f%-25.16f" % (50.0,50.0,50.0)

    model.grofile = gro_string

