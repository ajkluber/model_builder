import numpy as np

import mdtraj as md

from mdtraj.core.topology import Topology
from mdtraj.core.element import get_by_symbol

class CalphaMapping(object):
    """Calpha representation mapping"""
    def __init__(self, topology):
        self._source_topology = topology.copy()

        # Build new topology
        newTopology = Topology()
        new_atm_idx = 0 
        prev_ca = None
        ca_idxs = []
        for chain in topology._chains:
            newChain = newTopology.add_chain()
            for residue in chain._residues:
                resSeq = getattr(residue, 'resSeq', None) or residue.index
                newResidue = newTopology.add_residue(residue.name, newChain,
                                                     resSeq, residue.segment_id)
                # map CA
                new_ca = newTopology.add_atom('CA', get_by_symbol('C'), 
                                    newResidue, serial=new_atm_idx)

                ca_idxs.append([[ atm.index for atm in residue.atoms if \
                            (atm.name == "CA") ][0], new_atm_idx ])
                if prev_ca is None:
                    prev_ca = new_ca
                else:
                    newTopology.add_bond(prev_ca, new_ca)
                    prev_ca = new_ca
                new_atm_idx += 1

        self._ca_idxs = np.array(ca_idxs)
        self.topology = newTopology

    def map_traj(self, traj):
        """Create new trajectory"""
        ca_xyz = traj.xyz[:,self._ca_idxs[:,0],:]
        return md.Trajectory(ca_xyz, self.topology)

class CalphaCbetaMapping(object):
    """Calpha Cbeta center-of-mass representation mapping"""

    def __init__(self, topology):
        self._source_topology = topology.copy()

        newTopology = Topology()
        new_atm_idx = 0 
        prev_ca = None
        ca_idxs = []
        self._sidechain_idxs = []
        self._sidechain_mass = []
        for chain in topology._chains:
            newChain = newTopology.add_chain()
            for residue in chain._residues:
                resSeq = getattr(residue, 'resSeq', None) or residue.index
                newResidue = newTopology.add_residue(residue.name, newChain,
                                                     resSeq, residue.segment_id)
                # map CA
                new_ca = newTopology.add_atom('CA', get_by_symbol('C'), 
                                    newResidue, serial=new_atm_idx)
                if prev_ca is None:
                    prev_ca = new_ca
                else:
                    newTopology.add_bond(prev_ca, new_ca)
                    prev_ca = new_ca

                ca_idxs.append([[ atm.index for atm in residue.atoms if \
                            (atm.name == "CA") ][0], new_atm_idx ])
                new_atm_idx += 1

                if residue.name == 'GLY':
                    self._sidechain_idxs.append([])
                    self._sidechain_mass.append([])
                else:
                    # map CB
                    new_cb = newTopology.add_atom('CB', get_by_symbol('C'), 
                                        newResidue, serial=new_atm_idx)

                    newTopology.add_bond(new_cb, new_ca)

                    self._sidechain_idxs.append([[ atm.index for atm in residue.atoms if \
                                (atm.is_sidechain) and (atm.element.symbol != "H") ], new_atm_idx ])
                    self._sidechain_mass.append(np.array([ atm.element.mass for atm in residue.atoms if \
                                (atm.is_sidechain) and (atm.element.symbol != "H") ]))
                    new_atm_idx += 1

        self._ca_idxs = np.array(ca_idxs)
        self.topology = newTopology

    def map_traj(self, traj):
        """Return new Trajectory object with cacb topology and xyz"""
        cacb_xyz = np.zeros((traj.n_frames, self.topology.n_atoms, 3))

        for res in self.topology.residues:
            cacb_xyz[:,self._ca_idxs[res.index,1],:] = \
                    traj.xyz[:,self._ca_idxs[res.index,0],:]
            # Map sidechain atoms to their center of mass
            if res.name != 'GLY': 
                old_idxs = self._sidechain_idxs[res.index][0]
                new_idx = self._sidechain_idxs[res.index][1]
                sc_mass = self._sidechain_mass[res.index]
                tot_mass = np.sum(sc_mass)

                res_frms = traj.xyz[:,old_idxs,:]
                sc_com_xyz = np.array(map(lambda frm: \
                        np.sum(frm.T*sc_mass/tot_mass, axis=1), res_frms))

                cacb_xyz[:,new_idx,:] = sc_com_xyz

        return md.Trajectory(cacb_xyz, self.topology)

    
def write_bonds_tcl(topology, tcl_out="bonds.tcl"):
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
    for atm1, atm2 in topology.bonds:  
        tclstring += bondstring(molid,atm1.index, atm2.index) + "\n"

    with open(tcl_out, 'w') as fout:
        fout.write(tclstring)

if __name__ == "__main__":

    #name = "1JDP"
    name = "2KJV"
    traj = md.load(name+'.pdb')

    # test CA
    mapping = CalphaMapping(traj.top)
    ca_traj = mapping.map_traj(traj)
    ca_traj[0].save_pdb('ca_{}.pdb'.format(name))
    ca_traj.save_xtc('ca_{}.xtc'.format(name))

    # test CACB
    #mapping = CalphaCbetaMapping(traj.top)
    #cacb_traj = mapping.map_traj(traj)
    #cacb_traj[0].save_pdb('cacb_{}.pdb'.format(name))
    #cacb_traj.save_xtc('cacb_{}.xtc'.format(name))

    write_bonds_tcl(mapping.topology, tcl_out="{}_cabonds.tcl".format(name))

