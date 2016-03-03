from numpy import ndarray

from mdtraj import Topology, Trajectory

def write_bonds_tcl(bond_idxs, outfile="bonds.tcl"):

    bond_idxs = check_bond_idxs(bond_idxs)
    molid = 0
    bondstring = lambda molid, idx1, idx2: \
'''set sel [atomselect {0} "index {1} {2}"]
lassign [$sel getbonds] bond1 bond2
set id [lsearch -exact $bond1 {2}]
if {{ $id == -1 }} {{
lappend bond1 {2}
}}
set id [lsearch -exact $bond2 {1}]
if {{ $id == -1 }} {{
lappend bond2 {1}
}}
$sel setbonds [list $bond1 $bond2]
$sel delete'''.format(molid, idx1, idx2)

    tclstring = ''
    for idx1, idx2 in topology.bonds:  
        tclstring += bondstring(molid, idx1, idx2) + "\n"

    with open(outfile, 'w') as fout:
        fout.write(tclstring)

def write_bonds_conect(bond_idxs, outfile="conect.pdb"):
    
    bond_idxs = check_bond_idxs(bond_idxs)
 
    conectstring = ''
    for idx1, idx2 in bond_idxs:
        conectstring += "CONECT{:5d}{:5d}\n".format(idx1 + 1, idx2 + 1)

    with open(outfile, 'w') as fout:
        fout.write(conectstring)

def check_bond_idxs(thing):
    if type(thing) == Topology:
        bond_idxs = np.array([ [atm1.index, atm2.index] for atm1, atm2 in thing.bonds ])
    elif type(thing) == Trajectory:
        bond_idxs = np.array([ [atm1.index, atm2.index] for atm1, atm2 in thing.top.bonds ])
    elif type(thing) == ndarray:
        if topology.shape[1] != 2:
            raise IOError("array must be size (n_bonds, 2). Inputted: {}".format(topology.shape))
        else:
            bond_idxs = thing
         
    return bond_idxs

if __name__ == "__main__":
    import mdtraj as md
    traj = md.load("1SHG.pdb")
    write_bonds_tcl(traj.top)
    write_bonds_conect(traj.top)
