import numpy as np



############################################################################
# Calpha representation
#############################################################################
def _set_bonded_interactions(self):
    """ Extract info from the Native.pdb for making index and top file """
    ## Grab coordinates from the pdb file.
    self.cleanpdb = pdb_parser.get_clean_CA(self.pdb)
    self.cleanpdb_full = pdb_parser.get_clean_full(self.pdb)
    self.cleanpdb_full_noH = pdb_parser.get_clean_full_noH(self.pdb)
    coords, indices, atoms, residues = pdb_parser.get_coords_atoms_residues(self.cleanpdb)

    self.atom_indices = indices
    self.atom_types = atoms
    self.atom_residues = residues
    self.atom_coords = coords

    self.n_residues = len(residues)
    self.n_atoms = len(atoms)

    ## Set bonded force field terms quantities.
    self.bond_indices = [[indices[i],indices[i+1]] for i in range(self.n_atoms-1)]
    self.angle_indices = [[indices[i],indices[i+1],indices[i+2]] for i in range(self.n_atoms-2)]
    self.dihedral_indices = [[indices[i],indices[i+1],indices[i+2],indices[i+3]] for i in range(self.n_atoms-3)]

    self.bond_min = [ bond.distance(coords,i_idx-1,j_idx-1) for i_idx,j_idx in self.bond_indices ]
    self.angle_min = [ bond.angle(coords,i_idx-1,j_idx-1,k_idx-1) for i_idx,j_idx,k_idx in self.angle_indices ]
    self.dihedral_min = [ bond.dihedral(coords,i_idx-1,j_idx-1,k_idx-1,l_idx-1) for i_idx,j_idx,k_idx,l_idx in self.dihedral_indices ]

    if not hasattr(self,"bond_strengths"):
        self.bond_strengths = [ self.backbone_param_vals["Kb"] for i in range(len(self.bond_min)) ]
    if not hasattr(self,"angle_strengths"):
        self.angle_strengths = [ self.backbone_param_vals["Ka"] for i in range(len(self.angle_min)) ]
    if not hasattr(self,"dihedral_strengths"):
        self.dihedral_strengths = [ self.backbone_param_vals["Kd"] for i in range(len(self.dihedral_min)) ]




############################################################################
# Calpha Cbeta representation
#############################################################################
def _set_bonded_interactions(self):
    """ Extract info from the Native.pdb for making index and top file """
    ## Grab coordinates from the pdb file.
    self.cleanpdb = pdb_parser.get_clean_CA_CB(self.pdb)
    self.cleanpdb_full = pdb_parser.get_clean_full(self.pdb)
    self.cleanpdb_full_noH = pdb_parser.get_clean_full_noH(self.pdb)
    coords, indices, atoms, residues = pdb_parser.get_coords_atoms_residues(self.cleanpdb)

    self.atom_indices = indices
    self.atom_types = atoms
    self.atom_residues = residues
    self.atom_coords = coords

    self.n_residues = len(residues)
    self.n_atoms = len(atoms)

    ## Set bonded force field terms quantities.
    self.bond_indices = []
    for i in range(self.n_atoms-1):
        self.bond_indices.append([indices[2*i],indices[2*(i+1)]])
        self.bond_indices.append([indices[2*(i+1)],indices[2*(i+1)+1]])

    self.angle_indices = []
    for i in range(self.n_atoms-2):
        self.angle_indices.append([indices[i],indices[i+1],indices[i+2]])
        self.angle_indices.append([indices[i],indices[i+1],indices[i+2]])


    self.angle_indices = [[indices[i],indices[i+1],indices[i+2]] for i in range(self.n_atoms-2)]
    self.dihedral_indices = [[indices[i],indices[i+1],indices[i+2],indices[i+3]] for i in range(self.n_atoms-3)]

    self.bond_min = [ bond.distance(coords,i_idx-1,j_idx-1) for i_idx,j_idx in self.bond_indices ]
    self.angle_min = [ bond.angle(coords,i_idx-1,j_idx-1,k_idx-1) for i_idx,j_idx,k_idx in self.angle_indices ]
    self.dihedral_min = [ bond.dihedral(coords,i_idx-1,j_idx-1,k_idx-1,l_idx-1) for i_idx,j_idx,k_idx,l_idx in self.dihedral_indices ]

    if not hasattr(self,"bond_strengths"):
        self.bond_strengths = [ self.backbone_param_vals["Kb"] for i in range(len(self.bond_min)) ]
    if not hasattr(self,"angle_strengths"):
        self.angle_strengths = [ self.backbone_param_vals["Ka"] for i in range(len(self.angle_min)) ]
    if not hasattr(self,"dihedral_strengths"):
        self.dihedral_strengths = [ self.backbone_param_vals["Kd"] for i in range(len(self.dihedral_min)) ]
