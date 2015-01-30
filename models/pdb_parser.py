"""Clean and parse pdb files.

Find the full PDB format specification at:
http://www.wwpdb.org/documentation/format33/sect9.html#ATOM

PDB fixed-width column format is given by:
ATOM     44  C   ALA A  11      12.266  21.667  20.517  1.00 28.80           C  


To Do:
- Handle pdbs with multiple chains/molecules.
- Extend to parse other useful info from PDB files.
    - Calculate center of mass of sidechain atoms for use with
      Calpha-Cbeta model.

"""

import numpy as np


global atom_mass
atom_mass = {"H":1.00,"C":12.01,"N":14.00,"O":16.00,"S":32.07}

global backbone_atoms
backbone_atoms = ["N","CA","C","O","OXT"] 

#############################################################################
# Helper functions to clean the junk out of PDB files.
#############################################################################

def get_clean_CA(pdbname):
    """ Gets first chain from pdb. Keeps only CA atoms."""
    first = 0
    atomid = 1
    cleanpdb = ''
    for line in open(pdbname,'r'):
        line = line.rstrip("\n")
        if line[:3] in ['TER','END']:
            break
        else:
            # Keep only ATOM lines.
            if line[:4] == 'ATOM':
                if line[13:16].strip() == "CA":
                    if first == 0:
                        # Ignore alternative sidechain conformations.
                        if line[16] in ["A"," "]:
                            newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid,line[12:16],line[17:20],1,line[26:55])
                            atomid += 1
                            first = 1
                            first_index = int(line[22:26]) - 1
                            cleanpdb += newline
                    else:
                        # Ignore alternative sidechain conformations.
                        if line[16] in ["A"," "]:
                            newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid,line[12:16],line[17:20],int(line[22:26])-first_index,line[26:55])
                            atomid += 1
                            cleanpdb += newline
     
    cleanpdb += 'END\n'
    return cleanpdb

def get_clean_CA_CB(pdbname):
    """ Gets first chain from pdb. Keeps only CA CB atoms."""
    first = 0
    atomid = 1
    cleanpdb = ''
    for line in open(pdbname,'r'):
        line = line.rstrip("\n")
        if line[:3] in ['TER','END']:
            break
        else:
            # Keep only ATOM lines.
            if line[:4] == 'ATOM':
                if line[13:16].strip() in ["CA","CB"]:
                    if first == 0:
                        # Ignore alternative sidechain conformations.
                        if line[16] in ["A"," "]:
                            newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid,line[12:16],line[17:20],1,line[26:55])
                            atomid += 1
                            first = 1
                            first_index = int(line[22:26]) - 1
                            cleanpdb += newline
                    else:
                        # Ignore alternative sidechain conformations.
                        if line[16] in ["A"," "]:
                            newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid,line[12:16],line[17:20],int(line[22:26])-first_index,line[26:55])
                            atomid += 1
                            cleanpdb += newline
     
    cleanpdb += 'END\n'
    return cleanpdb

def get_clean_full(pdbname):
    """ Gets first chain from pdb. Keeps all atoms."""
    first = 0
    atomid = 1
    cleanpdb = ''
    for line in open(pdbname,'r'):
        line = line.rstrip("\n")
        if line[:3] in ['TER','END']:
            break
        else:
            # Keep only ATOM lines.
            if line[:4] == 'ATOM':
                if first == 0:
                    # Ignore alternative sidechain conformations.
                    if (line[16] in ["A"," "]) and (line[13] not in ["E","D"]):
                        newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                (atomid,line[12:16],line[17:20],1,line[26:55])
                        atomid += 1
                        first = 1
                        first_index = int(line[22:26]) - 1
                        cleanpdb += newline
                else:
                    # Ignore alternative sidechain conformations.
                    if (line[16] in ["A"," "]) and line[13] not in ["E","D"]:
                        newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                (atomid,line[12:16],line[17:20],int(line[22:26])-first_index,line[26:55])
                        atomid += 1
                        cleanpdb += newline
     
    cleanpdb += 'END\n'
    return cleanpdb

def get_clean_full_noH(pdbname):
    """ Gets first chain from pdb. Keeps all atoms except Hydrogen."""
    first = 0
    atomid = 1
    cleanpdb = ''
    for line in open(pdbname,'r'):
        line = line.rstrip("\n")
        if line[:3] in ['TER','END']:
            break
        else:
            # Keep only ATOM lines.
            if line[:4] == 'ATOM':
                if first == 0:
                    # Ignore alternative sidechain conformations.
                    if (line[16] in ["A"," "]) and (line[13] not in ["E","D"]):
                        # strip Hydrogens
                        if not line[12:16].strip().startswith("H"):
                            newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid,line[12:16],line[17:20],1,line[26:55])
                            atomid += 1
                            first = 1
                            first_index = int(line[22:26]) - 1
                            cleanpdb += newline
                else:
                    # Ignore alternative sidechain conformations.
                    if (line[16] in ["A"," "]) and line[13] not in ["E","D"]:
                        if not line[12:16].strip().startswith("H"):
                            newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid,line[12:16],line[17:20],int(line[22:26])-first_index,line[26:55])
                            atomid += 1
                            cleanpdb += newline
     
    cleanpdb += 'END\n'
    return cleanpdb

def get_coords_atoms_residues(pdb):
    """Parse lines of a pdb string.

    Parameters
    ----------
    pdb : str
        String in PDB file format that has contents of one protein chain. 
        Assumed to be one chain, numbering starting at 1. (e.g. as returned
        by get_clean_full_noH.
    pairs : array 
        Array of shape (n_pairs,2) that contains the atom indices to calculate
        the distance between.

    Returns
    -------
    atom_coords : array
        Array of atomic coordinates. Shape (n_atoms,3)
    atom_indices : array
        List of atom indices.
    atom_types : array
        List of atom types (C,CA,N,O,etc.).
    residue_indices : array
        List of residue indices.
    residue_types : array
        List of residue types (GLY,THR,etc.).
    residue_types_unique : array
        Nonredundant list of residue types (GLY,THR,etc.).

    See Also
    --------
    pdb_parser.get_clean_full_noH
    """
    atm_indxs = []
    atm_types = []
    atm_coords = []
    res_types = []
    res_types_unique = []
    res_indxs = []
    pdblines = pdb.split("\n")
    res_idx = 0
    for line in pdblines:
        if (line.startswith("END")) or (line.startswith("TER")):
            break
        else:
            atm_indxs.append(int(line[6:13]))
            atm_types.append(line[11:16].strip())
            atm_coords.append([float(line[31:39]),float(line[39:47]),float(line[47:55])]) 
            res_indxs.append(int(line[23:26]))
            res_types.append(line[17:20])
            if int(line[23:26]) == (res_idx + 1):
                res_types_unique.append(line[17:20])
                res_idx += 1

    # Coordinates in pdb files are Angstroms. Convert to nanometers.
    atm_coords = np.array(atm_coords)/10.
    atm_indxs = np.array(atm_indxs)
    atm_types = np.array(atm_types)
    res_indxs = np.array(res_indxs)
    res_types = np.array(res_types)

    return atm_coords,atm_indxs,atm_types,res_indxs,res_types,res_types_unique

def get_pairwise_distances(pdb,pairs):
    """Gets distances between pairs in pdb.

    Parameters
    ----------
    pdb : str
        String in PDB file format that has contents of one protein chain. 
        Assumed to be one chain, numbering starting at 1. (e.g. as returned
        by get_clean_full_noH.
    pairs : array 
        Array of shape (n_pairs,2) that contains the atom indices to calculate
        the distance between.

    Returns
    -------
    distances : array
        Array of distances between the pairs. Calculated in (nm).

    See Also
    --------
    pdb_parser.get_clean_full_noH
    """

    coords = get_coords_atoms_residues(pdb)[0]

    pairwise_distances = np.zeros(len(pairs),float)
    for i in range(len(pairs)):
        i_idx = pairs[i][0]
        j_idx = pairs[i][1]
        pairwise_distances[i] = np.linalg.norm(coords[i_idx-1] - coords[j_idx-1])
    return pairwise_distances

def calc_center_of_mass(atoms,coords):
    """Calculate center of mass of set of atoms"""
    com_xyz = np.zeros(3,float)
    n_atoms = len(atoms)
    total_mass = 0.
    for i in range(n_atoms):
        total_mass += atom_mass[atoms[i]]
        com_xyz += coords[i]*atom_mass[atoms[i]]
    com_xyz /= total_mass
    return com_xyz

def get_clean_CA_center_of_mass_CB(pdbname):
    """Get CA & CB of chain with CB at sidechain center of mass.

    Parameters
    ----------
    pdb : str
        String in PDB file format that has contents of one protein chain. 
        Assumed to be one chain, numbering starting at 1. (e.g. as returned
        by get_clean_full_noH.

    Returns
    -------
    cacab_pdb : str
        String in PDB format with only CA and CB atoms. With CB atoms placed
        at the sidechain center of mass.

    See Also
    --------
    pdb_parser.get_clean_full_noH
    """


    res_indx = 1
    atm_indx = 1
    coords = []
    atoms = []
    cacb_string = ""
    # Parse each line in the PDB string.
    for line in open(pdbname,"r").readlines():
        if (line.startswith("END")) or (line.startswith("TER")):
            # Process final residue.
            if prev_res_name != "GLY":
                # Skip glycines b/c they have no sidechain.
                com_xyz = calc_center_of_mass(atoms,coords)
                newline = 'ATOM%7d  %-4s%3s A%4d    %8.3f%8.3f%8.3f\n' % \
                    (atm_indx,"CB",prev_res_name,res_indx,com_xyz[0],com_xyz[1],com_xyz[2])
                atm_indx += 1
                cacb_string += newline
            break
        else:
            atom_type = line[11:16].strip()
            xyz = np.array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
            res_num = int(line[23:26])
            if atom_type == "CA":
                # We keep CA atoms just as they are.
                newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                        (atm_indx,line[12:16],line[17:20],res_num,line[26:55])
                atm_indx += 1
                cacb_string += newline

            if res_num == (res_indx + 1):
                # If we have moved to the next residuec calculate center of
                # mass of sidechain for previous residue.
                res_indx += 1 
                coords = np.array(coords)
                if prev_res_name != "GLY":
                    # Skip glycines b/c they have no sidechain.
                    com_xyz = calc_center_of_mass(atoms,coords)
                    newline = 'ATOM%7d  %-4s%3s A%4d    %8.3f%8.3f%8.3f\n' % \
                        (atm_indx,"CB",prev_res_name,res_indx-1,com_xyz[0],com_xyz[1],com_xyz[2])
                    atm_indx += 1
                    cacb_string += newline
                # Start collecting the next residue's info.
                if not atom_type in backbone_atoms: 
                    coords = [xyz]
                    atoms = [atom_type]
                else:
                    coords = []
                    atoms = []
            else:
                # Collect sidechain atom type and coords.
                prev_res_name = line[17:20]
                if not atom_type in backbone_atoms: 
                    coords.append(xyz)
                    atoms.append(atom_type[0])

    cacb_string += "END"

    return cacb_string

def get_CACB_contacts_from_AA_contact_map(pdbname,all_atom_map):
    """Get contact map for C-alpha C-beta model 

    Description
    -----------
    Determine the contact pairs based on a heavy-atom cutoff. 

    Parameters
    ----------
    pdbname : str
        The filename of a pdb.
    all_atom_map : str
        Filename of file that contains all-atom contact pairs.

    Returns
    -------
    CA_CA_pairs : list
        List of all contact C-alpha pairs. 
    CB_CB_pairs : list
        List of all contact C-beta pairs. 
    CA_CB_pairs : list
        List of all contact C-alpha C-beta pairs. 
    """

    aa_pairs = np.loadtxt(all_atom_map,dtype=int) 
    if len(aa_pairs[0,:]) == 4:
        aa_pairs = aa_pairs[:,1::2]

    pdb = get_clean_full_noH(pdbname)
    atm_coords,atm_indxs,atm_types,res_indxs,res_types,res_types_unique = get_coords_atoms_residues(pdb)
    cacb = get_clean_CA_center_of_mass_CB(pdbname)
    cacb_atm_indxs,cacb_atm_types,cacb_res_indxs = get_coords_atoms_residues(cacb)[1:4]
    n_res = len(np.unique(np.array(res_indxs,copy=True)))
    pairs = [[],[],[]]
    for n in range(len(aa_pairs)):
        # Determine which atoms are in contact. Which residues they are.
        atm1 = aa_pairs[n,0]
        atm2 = aa_pairs[n,1]
        convert_all_atom_contact_to_CACB(pairs,atm1,atm2,atm_types,res_indxs,cacb_res_indxs,cacb_atm_types)

    return pairs[0],pairs[1],pairs[2]

def convert_all_atom_contact_to_CACB(pairs,atm1,atm2,atm_types,res_indxs,cacb_res_indxs,cacb_atm_types):
    """Takes all-atom indices and determines the corresponding CA and/or CB indices"""
    atm1_type = atm_types[atm1 - 1]
    atm2_type = atm_types[atm2 - 1]
    atm1_resindx = res_indxs[atm1 - 1]
    atm2_resindx = res_indxs[atm2 - 1]
    res1 = (cacb_res_indxs == atm1_resindx).astype(int)
    res2 = (cacb_res_indxs == atm2_resindx).astype(int)
    ca = (cacb_atm_types == "CA").astype(int)
    cb = (cacb_atm_types == "CB").astype(int)

    if (atm1_type in backbone_atoms) and (atm2_type in backbone_atoms):
        # Create CA-CA contact.
        caindx1 = cacb_atm_indxs[(res1*ca).astype(bool)][0]
        caindx2 = cacb_atm_indxs[(res2*ca).astype(bool)][0]
        if [caindx1,caindx2] not in pairs[0]:
            pairs[0].append([caindx1,caindx2])
    elif (atm1_type in backbone_atoms) and (atm2_type not in backbone_atoms):
        # Create CB-CA contact
        caindx1 = cacb_atm_indxs[(res1*ca).astype(bool)][0]
        cbindx2 = cacb_atm_indxs[(res2*cb).astype(bool)][0]
        if [caindx1,cbindx2] not in pairs[1]:
            pairs[1].append([caindx1,cbindx2])
    elif (atm1_type not in backbone_atoms) and (atm2_type in backbone_atoms):
        # Create CA-CB contact
        cbindx1 = cacb_atm_indxs[(res1*cb).astype(bool)][0]
        caindx2 = cacb_atm_indxs[(res2*ca).astype(bool)][0]
        if [cbindx1,caindx2] not in pairs[1]:
            pairs[1].append([cbindx1,caindx2])
    else:
        # Create CB-CB contact
        cbindx1 = cacb_atm_indxs[(res1*cb).astype(bool)][0]
        cbindx2 = cacb_atm_indxs[(res2*cb).astype(bool)][0]
        if [cbindx1,cbindx2] not in pairs[2]:
            pairs[2].append([cbindx1,cbindx2])

def get_CACB_contacts_cutoff(pdbname,cutoff=0.45):
    """Get contact map for C-alpha C-beta model 

    Description
    -----------
    Determine the contact pairs based on a heavy-atom cutoff. 

    Parameters
    ----------
    pdbname : str
        The filename of a pdb.

    Returns
    -------
    pairs : list
        List of all contact pairs. 
    pairs_CA : list
        List of all contact C-alpha pairs. 
    pairs_CB : list
        List of all contact C-beta pairs. 
    pairs_CB_CA : list
        List of all contact C-alpha C-beta pairs. 
    """

    pdb = get_clean_full_noH(pdbname)
    atm_coords,atm_indxs,atm_types,res_indxs,res_types,res_types_unique = get_coords_atoms_residues(pdb)
    cacb = get_clean_CA_center_of_mass_CB(pdbname)
    cacb_info = get_coords_atoms_residues(cacb)
    n_res = len(np.unique(np.array(res_indxs,copy=True)))
    pairs = []
    pairs_CA = []
    pairs_CB = []
    pairs_CB_CA = []
    for i in range(1,n_res+1):
        res1info = (atm_coords[res_indxs == i],atm_types[res_indxs == i],atm_indxs[res_indxs == i])
        # Loop over residue pairs that are separated >= 4 in sequence
        for j in range(i+4,n_res+1):
            res2info = (atm_coords[res_indxs == j],atm_types[res_indxs == j],atm_indxs[res_indxs == j])
            crds2 = atm_coords[res_indxs == j]
            atms2 = atm_types[res_indxs == j]
            inds2 = atm_indxs[res_indxs == j]
            determine_contact(pairs,pairs_CA,pairs_CB,pairs_CB_CA,i,j,res1info,res2info,cacb_info,cutoff=cutoff) 
    return pairs,pairs_CA,pairs_CB,pairs_CB_CA

def determine_contact(pairs,pairs_CA,pairs_CB,pairs_CB_CA,i,j,res1info,res2info,cacb_info,cutoff=0.45):
    """If two atoms are within a cutoff determine what type of contact they're in."""
    # Unpack some inputs
    cacb_atm_indxs = cacb_info[1]
    cacb_atm_types = cacb_info[2]
    cacb_res_indxs = cacb_info[3]

    crds1 = res1info[0] ; crds2 = res2info[0]
    atms1 = res1info[1] ; atms2 = res2info[1]
    inds1 = res1info[2] ; inds2 = res2info[2]
    
    CB_CB = 0 ; CA_CA = 0
    CA_CB = 0 ; CB_CA = 0
    
    # If two residues have an atom within cutoff distance of
    # each other, assign one of these to the residue pair:
    for n in range(len(crds1)): 
        for m in range(len(crds2)): 
            dist = np.linalg.norm(crds1[n] - crds2[m])
            if dist <= cutoff:
                if (atms1[n] in backbone_atoms) and (atms2[m] in backbone_atoms):
                    # Backbone-backbone contact
                    CA_CA = 1
                elif (atms1[n] in backbone_atoms) and (atms2[m] not in backbone_atoms):
                    # Backbone-sidechain contact
                    CA_CB = 1
                elif (atms1[n] not in backbone_atoms) and (atms2[m] in backbone_atoms):
                    # Backbone-sidechain contact
                    CB_CA = 1
                elif (atms1[n] not in backbone_atoms) and (atms2[m] not in backbone_atoms):
                    # Sidechain-sidechain contact 
                    CB_CB = 1

    # If in contact, add the corresponding CA,CB atom indices. Also collect
    # whether these are sidechain-sidechain, mainchain-mainchain, or
    # mainchain-sidechain.
    if (CB_CB > 0) or (CB_CA > 0) or (CA_CB > 0) or (CA_CA > 0):
        res1 = (cacb_res_indxs == i).astype(int)
        res2 = (cacb_res_indxs == j).astype(int)
        ca = (cacb_atm_types == "CA").astype(int)
        cb = (cacb_atm_types == "CB").astype(int)
        if CB_CB > 0:
            if [i,j] not in pairs_CB:
                pairs_CB.append([i,j])
            cbindx1 = cacb_atm_indxs[(res1*cb).astype(bool)][0]
            cbindx2 = cacb_atm_indxs[(res2*cb).astype(bool)][0]
            if [cbindx1,cbindx2] not in pairs:
                pairs.append([cbindx1,cbindx2])

        if CA_CA > 0:
            if [i,j] not in pairs_CA:
                pairs_CA.append([i,j])
            caindx1 = cacb_atm_indxs[(res1*ca).astype(bool)][0]
            caindx2 = cacb_atm_indxs[(res2*ca).astype(bool)][0]
            if [caindx1,caindx2] not in pairs:
                pairs.append([caindx1,caindx2])

        if CB_CA > 0:
            if [i,j] not in pairs_CB_CA:
                pairs_CB_CA.append([i,j])
            cbindx1 = cacb_atm_indxs[(res1*cb).astype(bool)][0]
            caindx2 = cacb_atm_indxs[(res2*ca).astype(bool)][0]
            if [cbindx1,caindx2] not in pairs:
                pairs.append([cbindx1,caindx2])

        if CA_CB > 0:
            if [i,j] not in pairs_CB_CA:
                pairs_CB_CA.append([i,j])
            caindx1 = cacb_atm_indxs[(res1*ca).astype(bool)][0]
            cbindx2 = cacb_atm_indxs[(res2*cb).astype(bool)][0]
            if [caindx1,cbindx2] not in pairs:
                pairs.append([caindx1,cbindx2])

