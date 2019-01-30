from __future__ import absolute_import
import numpy as np

import mdtraj as md
from mdtraj.core.topology import Topology

import model_builder.models.mappings.util as util
import model_builder.models.mappings.atom_types as atom_types

class HeavyAtomMapping(object):

    def __init__(self, topology):
        r"""Calpha representation mapping

        Maps an all-atom representation to the heavy (non-hydrogen) atoms.

        Parameters
        ----------
        topology : mdtraj.Topology object

        """
        newTopology = md.Topology()     
        
        atom_mapping = {}

        atm_idx = 0
        res_idx = 1
        heavy_atom_idxs = []
        for chain in topology.chains:
            newChain = newTopology.add_chain()
            for residue in chain.residues:
                newResidue = newTopology.add_residue(residue.name, newChain, res_idx)
                for atom in residue.atoms:  
                    if atom.element.symbol not in ["H", "D"]:
                        new_atom = newTopology.add_atom(atom.name, 
                                            md.core.element.get_by_symbol(atom.element.symbol),
                                            newResidue, serial=atm_idx)
                        atom_mapping[atom] = new_atom
                        heavy_atom_idxs.append([atom.index, new_atom.index])
                        atm_idx += 1
                res_idx += 1
        
        # Add new bonds
        for atm1, atm2 in topology.bonds:
            if (atm1 in atom_mapping) and (atm2 in atom_mapping):
                new_atm1 = atom_mapping[atm1]
                new_atm2 = atom_mapping[atm2]
                newTopology.add_bond(new_atm1, new_atm2)

        self._heavy_atom_idxs = np.array(heavy_atom_idxs)
        self.topology = newTopology

    @property
    def top(self):
        return self.topology

    def map_traj(self, traj):
        """Create new trajectory"""
        hvy_xyz = traj.xyz[:,self._heavy_atom_idxs[:,0],:]
        return md.Trajectory(hvy_xyz, self.topology)
