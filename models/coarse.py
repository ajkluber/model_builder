import numpy as np
import argparse
import os
import shutil
import sys
import subprocess as sb

# I don't know yet if using arguments other than the standard, we'll see
def divide():
    top_file=open('smog.top')
    top_main=open('smog_main.top','w+')
    temp_top=open('temp.top','w+')
    for line in top_file:
        if len(line)>1 and line[0]<>';':
            if '[' in line and ']' in line:
                temp_top.close()
                section=line.strip().split(' ')[1]            
                temp_top=open('smog_%s.top'%section,'w+')
                top_main.write('[ %s ]\n'%section)
                top_main.write('#include "smog_%s.top"\n'%section)
                continue
            else:
                temp_top.write(line)
    
    temp_top.close()
    os.remove('temp.top')
    top_main.close()
    top_file.close()

def separate_dihedrals():
    original_dihedrals = open('smog_dihedrals.top')
    proper_dihedrals = open('smog_dihedrals_proper.top','w+')
    improper_dihedrals = open('smog_dihedrals_improper.top','w+')

    for line in original_dihedrals:
        if len(line.split())==7:
            proper_dihedrals.write('{0}'.format(line))
        elif len(line.split())==8:
            improper_dihedrals.write('{0}'.format(line))
        else:
            pass
    
    original_dihedrals.close()
    proper_dihedrals.close()
    improper_dihedrals.close()
    
    
def get_coarse_args():
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name',type=str, required=True, help='Name of protein')

    args = parser.parse_args()

    return args


# Read the atom list from the smog .ndx file
def read_ndx():
    index_file = open('smog.ndx','r').readlines()[1:]
    atom_list = np.zeros(len(index_file))

    for i in range(len(index_file)):
        atom_list[i] = index_file[i].split()[0]

    return atom_list

# Use the smog .top file to associate each atom with the residue it belongs to
def atoms_residues(atom_list):
    
    atom_res = np.zeros((len(atom_list),2))
    
    top_file = open('smog_atoms.top','r')
    
    for i in range(len(atom_list)):
        line = top_file.readline().split()
        atom_res[i][0] = int(line[0])
        atom_res[i][1] = int(line[2])

    return atom_res

# Read smog .contacts file
def read_contacts():

    contacts_file = np.loadtxt('smog.contacts')
    contacts_list = np.zeros((len(contacts_file),2))

    for i in range(len(contacts_file)):
        contacts_list[i][0] = int(contacts_file[i][1])
        contacts_list[i][1] = int(contacts_file[i][3])

    return contacts_list

# Translate contacts in terms of atom numbers to residue numbers
def determine_residue_contacts(atom_res, contacts_list):
    
    atoms = atom_res[:,0]
    residues = atom_res[:,1]

    residue_contacts = []

    for i in range(len(contacts_list)):
        a, b = contacts_list[i]
        a = int(a)
        b = int(b)
        c = np.where(atoms == a)
        d = np.where(atoms == b)
        residue_contacts.append((int(residues[c]), int(residues[d])))

    residue_contacts = np.array(residue_contacts)
    
    return residue_contacts

# Divide residues between long-range (|i-j|>8) and short-range
def trim_residues(residue_contacts):

    long_contacts = []
    short_contacts = []
    long_short_list = np.zeros(len(residue_contacts),dtype=int)

    for i in range(len(residue_contacts)):
        a, b = residue_contacts[i]

        if abs(a-b)>8:
            long_contacts.append((a,b))
            long_short_list[i] = 1
        else:
            short_contacts.append((a,b))
            long_short_list[i] = 0
        
    accumulated_pairs = []
    repeat_long_list = np.zeros(len(residue_contacts), dtype=int)

    for i in range(len(residue_contacts)):
        if long_short_list[i]==0:
            pass
        else:
            a,b = residue_contacts[i]
            if ((a,b) or (b,a)) not in accumulated_pairs:
                accumulated_pairs.append((a,b))
            else:
                repeat_long_list[i]=1
    
    long_contacts = np.array(long_contacts)
    short_contacts = np.array(short_contacts)
    accumulated_pairs = np.array(accumulated_pairs)

    np.savetxt("long_residue_contacts.dat",accumulated_pairs, fmt='%5d')

    return long_contacts, short_contacts, long_short_list, accumulated_pairs, repeat_long_list

# Build a list of residue numbers with the corresponding C_beta atom number
def get_c_beta_list():
    
    atom_file = open('smog_atoms.top','r').readlines()
    gro_file = open('smog.gro','r').readlines()

    c_beta_list = []

    for i in range(len(atom_file)):
        line = atom_file[i].split()
        coords = gro_file[i+2].split()
        if line[4] == "CB":
            c_beta_list.append((int(line[2]),int(line[0]),float(coords[3]),float(coords[4]),float(coords[5])))
            
    
    c_beta_file = open('c_beta_list.dat','w')
    for line in c_beta_list:
        c_beta_file.write('{0:4d}    {1:4d}   {2:3.3f}    {3:3.3f}    {4:3.3f}\n'.format(line[0],line[1],line[2],line[3],line[4]))
    
    c_beta_file.close()

    return c_beta_list
    
def determine_c_beta_distances(atom_1, atom_2):
    
    a = np.array([atom_1[2],atom_1[3],atom_1[4]])
    b = np.array([atom_2[2],atom_2[3],atom_2[4]])

    distance = np.linalg.norm(a-b)

    return distance

    

def modify_pairs_exclusions(residue_contacts, long_short_list, repeat_long_list, c_beta_list):
    
    pairs_file = np.loadtxt('smog_pairs.top')
    exclusions_file = np.loadtxt('smog_exclusions.top')

    shutil.move('smog_pairs.top', 'smog_pairs_top_old')
    shutil.move('smog_exclusions.top','smog_exclusions_top_old')

    #Differentiate between long and short contact pairs

    short_pairs_file = open('smog_pairs_s.top','w')
#    fixed_long_pairs_file_1 = open('smog_pairs_f1.top','w')
#    fixed_long_pairs_file_2 = open('smog_pairs_f2.top','w')
    long_pairs_file = open('smog_pairs_long','w')
    long_pairs_top = open('smog_pairs_l.top','w')
    long_bonds_rep = open('smog_bonds_rep.top','w')
    new_exclusions_file = open('smog_exclusions.top','w')

    residue_list = np.array([row[0] for row in c_beta_list])
    print residue_list

    # Minimum of well for Gaussian contacts (sigma) is column 4
    bond_rep = 0
    # Multiplicity factor is the factor by which all long range coarse epsilons are calculated in the vanilla model so as to preserve the energy balance given by smog2

    all_long_bonds = np.sum(np.array(long_short_list))
    coarse_long_bonds = 0
    for i in range(len(long_short_list)):
        if repeat_long_list[i]==1:
            pass
        else:
            if long_short_list[i]==1:
                res_1 = residue_contacts[i,0]
                res_2 = residue_contacts[i,1]

                lookup_a = np.where(residue_list == res_1)[0]
                lookup_b = np.where(residue_list == res_2)[0]

                if lookup_a.size == 0:
                    pass
                elif lookup_b.size == 0:
                    pass
                else:
                    coarse_long_bonds+=1
                    
    multiplicity_factor = float(all_long_bonds)/float(coarse_long_bonds)
 
    for i in range(len(long_short_list)):

        if repeat_long_list[i]==1:
            pass
        else:

            pair = pairs_file[i]
            exc = exclusions_file[i]

            if long_short_list[i]==0:
                
                short_pairs_file.write('{0:4d}   {1:4d}    {2:2d}    {3:2.12e}   {4:2.12e}   {5:2.12e}   {6:2.12e}\n'.format(int(pair[0]),int(pair[1]),int(pair[2]),pair[3],pair[4],pair[5],pair[6]))
                new_exclusions_file.write('{0:4d}   {1:4d}\n'.format(int(exc[0]),int(exc[1])))

            elif long_short_list[i]==1:
                # 1) Correlate residues with their corresponding C_beta number and coordinates
                # 2) Calculate new equilibrium distance
                # 3) Write new pair and inclusion file lines
                res_1 = residue_contacts[i,0]
                res_2 = residue_contacts[i,1]
                
                lookup_a = np.where(residue_list == res_1)[0]
                lookup_b = np.where(residue_list == res_2)[0]

                if lookup_a.size == 0:
                    print 'Residue {0} is a GLY and therefore loses its contacts in a C_beta model'.format(res_1)
                    pass
                elif lookup_b.size == 0:
                    print 'Residue {0} is a GLY and therefore loses its contacts in a C_beta model'.format(res_2)
                    pass
                else:
                    a = int(lookup_a)
                    b = int(lookup_b)
    
                    atom_1 = c_beta_list[a]
                    atom_2 = c_beta_list[b]
                
                    new_r0 = determine_c_beta_distances(atom_1,atom_2)
                    new_sigma = new_r0/np.sqrt(50*np.log(2))
                    new_epsilon = pair[3]*multiplicity_factor
                    long_pairs_top.write('{0:4d}   {1:4d}    {2:2d}    {3:2.12e}   {4:2.12e}   {5:2.12e}\n'.format(int(atom_1[1]),int(atom_2[1]),int(5),new_epsilon,new_r0,new_sigma))
                    
                    #Eliminating the exclusions for the long range contacts. Will use the standard LJ12 repulsive wall
                   # new_exclusions_file.write('{0:4d}   {1:4d}\n'.format(int(atom_1[1]),int(atom_2[1])))

                    # This file will keep the original list of pairs.
                    long_pairs_file.write('{0:4d}   {1:4d}    {2:2d}    {3:2.12e}   {4:2.12e}   {5:2.12e}\n'.format(int(atom_1[1]),int(atom_2[1]),int(5),new_epsilon,new_r0,new_sigma))

                    # I will consider that the third term of the Lammert et al. model is negligible since the excluded radius is much smaller than the equilibrium contact radius for the long range interactions.

#                    fixed_long_pairs_file_1.write('{0:4d}   {1:4d}    {2:2d}    {3:2.12e}   {4:2.12e}   {5:2.12e}   {6:2.12e}\n'.format(int(atom_1[1]),int(atom_2[1]),int(6),new_epsilon,new_r0,new_sigma,pair[6]))

#                    fixed_long_pairs_file_2.write('{0:4d}   {1:4d}    {2:2d}    {3:2.12e}   {4:2.12e}   {5:2.12e}\n'.format(int(atom_1[1]),int(atom_2[1]),int(5),-pair[3],new_r0,new_sigma))

    short_pairs_file.close()
    long_pairs_file.close()
    long_pairs_top.close()
    long_bonds_rep.close()
    new_exclusions_file.close()
                                                  
    
    main_file = open('smog_main.top','r')
    new_main_file = open('smog_main_n.top','w')
    for line in main_file:
        if line.split()[1] == '"smog_pairs.top"':
            new_main_file.write('#include "smog_pairs_s.top"\n')
            new_main_file.write('#include "smog_pairs_l.top"\n')
#            new_main_file.write('#include "smog_pairs_f1.top"\n')
#            new_main_file.write('#include "smog_pairs_f2.top"\n')
                                                
        elif line.split()[1] == '"smog_bonds.top"':
            new_main_file.write('#include "smog_bonds.top"\n')
            new_main_file.write('#include "smog_bonds_rep.top"\n')
        else:
            new_main_file.write(line)

    new_main_file.close()
            
    os.remove('smog_main.top')
    shutil.move('smog_main_n.top','smog_main.top')
        
    
def generate_long_short_index():
    
    long_pairs = np.loadtxt('smog_pairs_long')
    short_pairs = np.loadtxt('smog_pairs_s.top')
    
    long_pairs = long_pairs[:,:2]
    short_pairs = short_pairs[:,:2]

    long_file = open('smog_long.ndx','w')
    long_file.write('[pairs]\n')
    short_file = open('smog_short.ndx','w')
    short_file.write('[pairs]\n')
    all_file = open('smog_long_and_short.ndx','w')
    all_file.write('[pairs]\n')

    for i in range(len(long_pairs)):
        long_file.write('{0:5d}   {1:5d}\n'.format(int(long_pairs[i,0]),int(long_pairs[i,1])))
        all_file.write('{0:5d}   {1:5d}\n'.format(int(long_pairs[i,0]),int(long_pairs[i,1])))
   
    for i in range(len(short_pairs)):
        short_file.write('{0:5d}   {1:5d}\n'.format(int(short_pairs[i,0]),int(short_pairs[i,1])))
        all_file.write('{0:5d}   {1:5d}\n'.format(int(short_pairs[i,0]),int(short_pairs[i,1])))
    
    long_file.close()
    short_file.close()
    all_file.close()

def generate_files(directory):
    
    pdb_line = 'editconf -f smog.gro -n smog.ndx -o smog.pdb'
    sb.call(pdb_line.split())

    divide()

    separate_dihedrals()
    
    atom_list = read_ndx()

    atom_res = atoms_residues(atom_list)
    
    contacts_list = read_contacts()
    
    residue_contacts = determine_residue_contacts(atom_res, contacts_list)

    long_contacts, short_contacts, long_short_list, accumulated_pairs, repeat_long_list = trim_residues(residue_contacts)
    
    c_beta_list = get_c_beta_list()
    
    modify_pairs_exclusions(residue_contacts, long_short_list, repeat_long_list, c_beta_list)

    generate_long_short_index()

if __name__ == "__main__":
    
    directory = get_coarse_args()

    generate_files(directory)
