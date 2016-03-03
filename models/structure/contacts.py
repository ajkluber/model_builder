import itertools
import numpy as np

def residue_contacts(traj, cutoff=0.5, exclude_neighbors=4):
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
    return contacts
