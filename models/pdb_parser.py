''' Parse and clean pdb file.

Find the full PDB format specification at:
http://www.wwpdb.org/documentation/format33/sect9.html#ATOM

PDB fixed-width column format is given by:
ATOM     44  C   ALA A  11      12.266  21.667  20.517  1.00 28.80           C  




To Do:
- Handle pdbs with multiple chains/molecules.
- Extend to parse other useful info from PDB files.
    - Calculate center of mass of sidechain atoms for use with
      Calpha-Cbeta model.

'''

import numpy as np

import bonded_potentials as bond

global atom_mass
atom_mass = {"C":12.01,"N":14.00,"O":16.00,"S":32.07}

def get_clean_CA(pdbname):
    ''' Gets first chain from pdb. Keeps only CA atoms.'''
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
    ''' Gets first chain from pdb. Keeps only CA CB atoms.'''
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
    ''' Gets first chain from pdb. Keeps all atoms.'''
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
    ''' Gets first chain from pdb. Keeps all atoms except Hydrogen.'''
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
    ''' Parse lines of a pdb string. Returns coordinates in nm. '''
    indices = []
    atoms = []
    residues = []
    coords = []
    pdblines = pdb.split("\n")
    res_indx = 0
    for line in pdblines:
        if line.startswith("END"):
            break
        else:
            indices.append(int(line[6:13]))
            atoms.append(line[11:16].strip())
            coords.append([float(line[31:39]),float(line[39:47]),float(line[47:55])]) 
            if (int(line[23:26]) == (res_indx + 1)):
                res_indx += 1 
                residues.append(line[17:20])

    ## Coordinates in pdb files are Angstroms. Convert to nanometers.
    coords = np.array(coords)/10.

    return coords,indices,atoms,residues



def get_pairwise_distances(pdb,pairs):
    ''' Calculate atomic distances between pairs in nm. '''

    coords,indices,atoms,residues = get_coords_atoms_residues(pdb)

    pairwise_distances = np.zeros(len(pairs),float)
    for i in range(len(pairs)):
        i_idx = pairs[i][0]
        j_idx = pairs[i][1]
        pairwise_distances[i] = bond.distance(coords,i_idx-1,j_idx-1)
    return pairwise_distances

def get_clean_CA_center_of_mass_CB(pdb):
    """ Return pdb with CACB atoms but with CB's at sidechain center of mass """

    backbone_atoms = ["N","CA","C","O"]

    pdblines = pdb.split("\n")
    res_indx = 1
    atm_indx = 1
    sidechain_coords = []
    sidechain_atoms = []
    cacb_string = ""
    for line in pdblines:
        if line.startswith("END"):
            if prev_res_name != "GLY":
                com_xyz = np.zeros(3,float)
                n_atoms = len(sidechain_atoms)
                total_mass = 0.
                for i in range(n_atoms):
                    total_mass += atom_mass[sidechain_atoms[i]]
                    com_xyz += sidechain_coords[i]*atom_mass[sidechain_atoms[i]]
                com_xyz /= total_mass

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
                newline = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                        (atm_indx,line[12:16],line[17:20],res_num,line[26:55])
                atm_indx += 1
                cacb_string += newline

            if res_num == (res_indx + 1):
                ## Calculate center of mass of sidechain for previous residue.
                res_indx += 1 
                sidechain_coords = np.array(sidechain_coords)
                if prev_res_name != "GLY":
                    com_xyz = np.zeros(3,float)
                    n_atoms = len(sidechain_atoms)
                    total_mass = 0.
                    for i in range(n_atoms):
                        total_mass += atom_mass[sidechain_atoms[i]]
                        com_xyz += sidechain_coords[i]*atom_mass[sidechain_atoms[i]]
                    com_xyz /= total_mass

                    newline = 'ATOM%7d  %-4s%3s A%4d    %8.3f%8.3f%8.3f\n' % \
                            (atm_indx,"CB",prev_res_name,res_indx-1,com_xyz[0],com_xyz[1],com_xyz[2])
                    atm_indx += 1
                    cacb_string += newline
                if not atom_type in backbone_atoms: 
                    sidechain_coords = [xyz]
                    sidechain_atoms = [atom_type]
                else:
                    sidechain_coords = []
                    sidechain_atoms = []
            else:
                prev_res_name = line[17:20]
                if not atom_type in backbone_atoms: 
                    sidechain_coords.append(xyz)
                    sidechain_atoms.append(atom_type[0])

    cacb_string += "END"

    return cacb_string
