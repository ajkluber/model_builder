import itertools
import numpy as np


def heavy_atom_cutoff(ires, jres, ref_xyz, contacts, cutoff):
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

def cb_cutoff(ires, jres, ref_xyz, contacts, cutoff):
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
        contacts.append((ires.index, jres.index))

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
                heavy_atom_cutoff(ires, jres, ref_xyz, contacts, cutoff)
            elif method == "cb_cutoff":
                cb_cutoff(ires, jres, ref_xyz, contacts, cutoff)
            else:
                raise ValueError("method arg should be heavy_atom_cutoff or cb_cutoff")

    return contacts
