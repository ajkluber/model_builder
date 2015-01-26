""" Parse and clean pdb file.

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
            ## Keep only ATOM lines.
            if line[:4] == 'ATOM':
                if line[13:16].strip() == "CA":
                    if first == 0:
                        if line[16] in ["A"," "]:
                            newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid,line[12:16],line[17:20],1,line[26:55])
                            atomid += 1
                            first = 1
                            first_index = int(line[22:26]) - 1
                            cleanpdb += newline
                    else:
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
            ## Keep only ATOM lines.
            if line[:4] == 'ATOM':
                if line[13:16].strip() in ["CA","CB"]:
                    if first == 0:
                        if line[16] in ["A"," "]:
                            newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid,line[12:16],line[17:20],1,line[26:55])
                            atomid += 1
                            first = 1
                            first_index = int(line[22:26]) - 1
                            cleanpdb += newline
                    else:
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
            ## Keep only ATOM lines.
            if line[:4] == 'ATOM':
                if first == 0:
                    if (line[16] in ["A"," "]) and (line[13] not in ["E","D"]):
                        newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                (atomid,line[12:16],line[17:20],1,line[26:55])
                        atomid += 1
                        first = 1
                        first_index = int(line[22:26]) - 1
                        cleanpdb += newline
                else:
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
            ## Keep only ATOM lines.
            if line[:4] == 'ATOM':
                if first == 0:
                    if (line[16] in ["A"," "]) and (line[13] not in ["E","D"]):
                        ## strip Hydrogens
                        if not line[12:16].strip().startswith("H"):
                            newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid,line[12:16],line[17:20],1,line[26:55])
                            atomid += 1
                            first = 1
                            first_index = int(line[22:26]) - 1
                            cleanpdb += newline
                else:
                    if (line[16] in ["A"," "]) and line[13] not in ["E","D"]:
                        if not line[12:16].strip().startswith("H"):
                            newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid,line[12:16],line[17:20],int(line[22:26])-first_index,line[26:55])
                            atomid += 1
                            cleanpdb += newline
     
    cleanpdb += 'END\n'
    return cleanpdb

def get_coords_atoms_residues(pdb):
    """ Parse lines of a pdb string. Returns coordinates in nm. """
    atm_indxs = []
    atm_types = []
    atm_coords = []
    res_types = []
    res_indxs = []
    pdblines = pdb.split("\n")
    res_indx = 0
    for line in pdblines:
        if line.startswith("END"):
            break
        else:
            atm_indxs.append(int(line[6:13]))
            atm_types.append(line[11:16].strip())
            atm_coords.append([float(line[31:39]),float(line[39:47]),float(line[47:55])]) 
            if (int(line[23:26]) == (res_indx + 1)):
                res_indxs.append(int(line[23:26]))
                res_types.append(line[17:20])
                res_indx += 1 

    ## Coordinates in pdb files are Angstroms. Convert to nanometers.
    atm_coords = np.array(atm_coords)/10.

    return atm_coords,atm_indxs,atm_types,res_indxs,res_types

def get_pairwise_distances(pdb,pairs):
    """Calculate atomic distances between pairs in nm."""

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

def get_clean_CA_center_of_mass_CB(pdb):
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

    backbone_atoms = ["N","CA","C","O"] 

    pdblines = pdb.split("\n")
    res_indx = 1
    atm_indx = 1
    coords = []
    atoms = []
    cacb_string = ""
    # Parse each line in the PDB string.
    for line in pdblines:
        if (line.startswith("END")) or (line.startswith("TER")):
            # Process final residue.
            if prev_res_name != "GLY":
                # Skip glycines b/c they have no sidechain.
                com_xyz = calc_center_of_mass(atoms,coords)
                newline = 'ATOM%7d  %-4s%3s A%4d    %8.3f%8.3f%8.3f\n' % \
                    (atm_indx,"CB",prev_res_name,res_indx-1,com_xyz[0],com_xyz[1],com_xyz[2])
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
