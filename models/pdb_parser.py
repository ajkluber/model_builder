''' Parse and clean pdb file.

Find the full PDB format specification at:
http://www.wwpdb.org/documentation/format33/sect9.html#ATOM

PDB fixed-width column format is given by:
ATOM     44  C   ALA A  11      12.266  21.667  20.517  1.00 28.80           C  




To Do:
- Handle pdbs with multiple chains/molecules.
- Extend to parse other useful info from PDB files.

'''

import numpy as np

import bonded_potentials as bond

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
