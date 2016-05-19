import itertools
import numpy as np


def _heavy_atom_cutoff(ires, jres, ref_xyz, contacts, cutoff):
    # Calculate the minimum distance between heavy atoms
    iheavy_xyz = ref_xyz[[ atm.index for atm in ires.atoms \
                    if atm.element.symbol != 'H' ],:]
    jheavy_xyz = ref_xyz[[ atm.index for atm in jres.atoms \
                    if atm.element.symbol != 'H' ],:]
    min_dist = np.min([ np.linalg.norm(n_xyz - k_xyz) \
                        for n_xyz in iheavy_xyz \
                        for k_xyz in jheavy_xyz ])
    if min_dist <= cutoff:
        contacts.append((ires.index, jres.index))
        # To Do: 
        # Shadow algorithm would require another loop right here
        # to check for intermediate residues

def _cb_cutoff(ires, jres, ref_xyz, contacts, cutoff):
    # Calculate the CB cutoff contact
    if ires.name == "GLY":
        cb1_xyz = ref_xyz[[ atom.index for atom in  ires.atoms if atom.name == "CA" ],:][0]
    else:
        cb1_xyz = ref_xyz[[ atom.index for atom in  ires.atoms if atom.name == "CB" ],:][0]
 
    if jres.name == "GLY":
        cb2_xyz = ref_xyz[[ atom.index for atom in  jres.atoms if atom.name == "CA" ],:][0]
    else:
        cb2_xyz = ref_xyz[[ atom.index for atom in  jres.atoms if atom.name == "CB" ],:][0]
 
    dist = np.linalg.norm(cb1_xyz - cb2_xyz)

    if dist <= cutoff:
        contacts.append((ires, jres))

def residue_contacts(traj, cutoff=0.5, exclude_neighbors=4, method="heavy_atom_cutoff"):
    """Return indices of residues that are in contact"""
    if traj.n_frames > 1:
        print "warning: only using first frame of as reference"
    ref_xyz = traj[0].xyz[0] 

    contacts = []
    for ires, jres in itertools.product(traj.top.residues, traj.top.residues):
        if (ires.chain.index == jres.chain.index) and \
            (abs(ires.index - jres.index) < exclude_neighbors):
            # Exclude neighbors on same chain
            continue
        elif jres.index > ires.index:
            if method == "heavy_atom_cutoff":
                _heavy_atom_cutoff(ires, jres, ref_xyz, contacts, cutoff)
            elif method == "cb_cutoff":
                _cb_cutoff(ires, jres, ref_xyz, contacts, cutoff)
            else:
                raise ValueError("method arg should be heavy_atom_cutoff or cb_cutoff")

    return contacts

def _standard_resname_and_charge(res_name):
    """Convert non-standard residue names and return charge"""
    std_resnames = {"ARN":"ARG", "ASH":"ASP", "GLH":"GLU", "LYN":"LYS", "HIE":"HIS", "HIP":"HIS"}
    std_charges = {"ARG":1, "ASP":-1, "GLU":-1, "LYS":1, "HIE":1, "HIP":2, "HIS":1}

    # determine charge from name
    if res_name in std_charges.keys():
        charge = std_charges[res_name]
    else:
        charge = 0

    # change name if needed.
    if res_name in std_resnames.keys():
        new_res_name = std_resnames[res_names] 
    else:
        new_res_name = res_name
    
    return new_res_name, charge


def trilaterate(v1, v2, v3, d1, d2, d3, left_handed=True):
    """Find position of fourth atom with three distances of three other atoms.

    Parameters
    ----------
    v1 : np.ndarray(3) 
        Coordinates of first atom. Used as origin.
    v2 : np.ndarray(3) 
        Coordinates of second atom. Aligned to x-axis.
    v3 : np.ndarray(3) 
        Coordinates of third atom. Placed in the xy-plane.
    d1 : float
        Distance to first atom.
    d2 : float
        Distance to second atom.
    d3 : float
        Distance to third atom.

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Trilateration
    """

    # define orthonormal basis vectors e_x, e_y, e_z
    v2_x = np.linalg.norm(v2 - v1)
    e_x = (v2 - v1)/v2_x

    v3_x = np.dot(e_x, v3 - v1)

    e_y = (v3 - v1 - v3_x*e_x)/np.linalg.norm(v3 - v1 - v3_x*e_x)

    v3_y = np.dot(e_y, v3 - v1)

    e_z = np.cross(e_x, e_y)

    # Should be unit norm.
    #print np.linalg.norm(e_x), np.linalg.norm(e_y), np.linalg.norm(e_z) 

    v4_x = (d1**2 - d2**2 + v2_x**2)/(2.*v2_x)
    v4_y = ((d1**2 - d3**2 + v3_x**2 + v3_y**2)/(2.*v3_y)) - (v3_x/v3_y)*v4_x
    v4_z = np.sqrt(d1**2 - v4_x**2 - v4_y**2)

    # two roots indicate left- and right-handed chiralities. 
    if left_handed: 
        # amino acid sidechains are left-handed.
        v4 = v1 + v4_x*e_x + v4_y*e_y - v4_z*e_z
    else:
        v4 = v1 + v4_x*e_x + v4_y*e_y + v4_z*e_z

    return v4
