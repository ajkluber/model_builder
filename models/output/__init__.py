

class GromacsSimulationFiles(object):


    def __init__(self, Hamiltonian, version=None):
        self.Hamiltonian
        pass

    # Check pair opts


def get_atoms_string(Hamiltonian):
    """ Generate the [ atoms ] string."""
    atoms_string = " [ atoms ]\n"
    atoms_string += " ;nr  type  resnr residue atom  cgnr charge  mass\n"
    for j in range(len(Hamiltonian.atm_indxs)):
        atmnum = Hamiltonian.atm_indxs[j]
        atmtype = Hamiltonian.atm_names[j] # changed atm_types to atm_names
        resnum = Hamiltonian.res_indxs[j]
        restype = Hamiltonian.res_types[j]
        atoms_string += " %5d%4s%8d%5s%4s%8d%8.3f%8.3f\n" % \
                    (atmnum,atmtype,resnum,restype,atmtype,atmnum,0.0,1.0)
    Hamiltonian.atoms_string = atoms_string

def get_bonds_string(Hamiltonian):
    """ Generate the [ bonds ] string."""
    bonds_string = " [ bonds ]\n"
    bonds_string += " ; ai aj func r0(nm) Kb\n"
    for k in range(len(Hamiltonian.bond_min)):
        i_idx = Hamiltonian.bond_indices[k][0]
        j_idx = Hamiltonian.bond_indices[k][1]
        dist = Hamiltonian.bond_min[k]
        kb = Hamiltonian.bond_strengths[k]
        bonds_string += "%6d %6d%2d%18.9e%18.9e\n" %  \
                      (i_idx,j_idx,1,dist,kb)
    Hamiltonian.bonds_string = bonds_string

def get_angles_string(Hamiltonian):
    """ Generate the [ angles ] string."""
    angles_string = " [ angles ]\n"
    angles_string += " ; ai  aj  ak  func  th0(deg)   Ka\n"
    for n in range(len(Hamiltonian.angle_min)):
        i_idx = Hamiltonian.angle_indices[n][0]
        j_idx = Hamiltonian.angle_indices[n][1]
        k_idx = Hamiltonian.angle_indices[n][2]
        theta = Hamiltonian.angle_min[n]
        ka = Hamiltonian.angle_strengths[n]
        angles_string += "%6d %6d %6d%2d%18.9e%18.9e\n" %  \
                      (i_idx,j_idx,k_idx,1,theta,ka)
    Hamiltonian.angles_string = angles_string

def get_dihedrals_string(Hamiltonian):
    """ Generate the [ dihedrals ] string."""
    dihedrals_string = " [ dihedrals ]\n"
    dihedrals_string += " ; ai  aj  ak al  func  phi0(deg)   Kd mult\n"
    dihedrals_ndx = '[ dihedrals ]\n'
    for n in range(len(Hamiltonian.dihedral_min)):
        i_idx = Hamiltonian.dihedral_indices[n][0]
        j_idx = Hamiltonian.dihedral_indices[n][1]
        k_idx = Hamiltonian.dihedral_indices[n][2]
        l_idx = Hamiltonian.dihedral_indices[n][3]
        dih_type = Hamiltonian.dihedral_type[n]
        phi = Hamiltonian.dihedral_min[n]
        kd = Hamiltonian.dihedral_strengths[n]
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
    Hamiltonian.dihedrals_string = dihedrals_string
    Hamiltonian.dihedrals_ndx = dihedrals_ndx
