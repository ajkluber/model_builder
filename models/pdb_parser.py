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

def clean(pdbname):
    ''' Filter pdb lines. Ignores everything before first ATOM line and after first END or TER '''
    first_full = 0
    atomid_full = 1
    cleanpdb_full = ''
    cleanpdb_full_noH = ''
    first_ca = 0
    atomid_ca = 1
    cleanpdb_ca = ''
    for line in open(pdbname,'r'):
        line = line.rstrip("\n")
        if line[:3] in ['TER','END']:
            break
        else:
            ## Keep only ATOM lines.
            if line[:4] == 'ATOM':
                if line[13:16].strip() == "CA":
                    if first_ca == 0:
                        if line[16] in ["A"," "]:
                            newline_ca = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid_ca,line[12:16],line[17:20],1,line[26:55])
                            atomid_ca += 1
                            first_ca = 1
                            first_index_ca = int(line[22:26]) - 1
                            cleanpdb_ca += newline_ca
                    else:
                        if line[16] in ["A"," "]:
                            newline_ca = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                    (atomid_ca,line[12:16],line[17:20],int(line[22:26])-first_index_ca,line[26:55])
                            atomid_ca += 1
                            cleanpdb_ca += newline_ca

                if first_full == 0:
                    if (line[16] in ["A"," "]) and (line[13] not in ["E","D"]):
                        newline_full = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                (atomid_full,line[12:16],line[17:20],1,line[26:55])
                        atomid_full += 1
                        first_full = 1
                        first_index_full = int(line[22:26]) - 1
                        cleanpdb_full += newline_full
                        ## strip Hydrogens
                        if not line[12:16].strip().startswith("H"):
                            cleanpdb_full_noH += newline_full
                else:
                    if (line[16] in ["A"," "]) and line[13] not in ["E","D"]:
                        newline_full = 'ATOM%7s %-5s%3s A%4d%s\n' % \
                                (atomid_full,line[12:16],line[17:20],int(line[22:26])-first_index_full,line[26:55])
                        atomid_full += 1
                        cleanpdb_full += newline_full
                        if not line[12:16].strip().startswith("H"):
                            cleanpdb_full_noH += newline_full
     
    cleanpdb_full += 'END\n'
    cleanpdb_full_noH += 'END\n'
    cleanpdb_ca += 'END\n'
    return cleanpdb_full,cleanpdb_full_noH,cleanpdb_ca

def get_coords_atoms_residues(pdb):
    ''' Parse lines of a pdb string'''
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
