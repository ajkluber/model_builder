import numpy as np

import mdtraj as md
from mdtraj.core.topology import Topology

import util
import atom_types

class CalphaMapping(object):
    r"""Calpha representation mapping"""

    def __init__(self, topology):
        r"""Calpha representation mapping

        Maps an all-atom representation to just the C-alpha's of the backbone.

        Holds default assignment of .

        Parameters
        ----------
        topology : mdtraj.Topology object

        """


        self._ref_topology = topology.copy()

        # Build new topology
        newTopology = Topology()
        prev_ca = None
        ca_idxs = []
        atm_idx = 0 
        for chain in topology._chains:
            newChain = newTopology.add_chain()
            for residue in chain._residues:
                resSeq = getattr(residue, 'resSeq', None) or residue.index
                newResidue = newTopology.add_residue(residue.name, newChain, resSeq)
                # map CA
                new_ca = newTopology.add_atom('CA', md.core.element.get_by_symbol('C'), 
                                    newResidue, serial=atm_idx)

                ca_idxs.append([[ atm.index for atm in residue.atoms if \
                            (atm.name == "CA") ][0], atm_idx ])
                if prev_ca is None:
                    prev_ca = new_ca
                else:
                    if prev_ca.residue.chain.index == new_ca.residue.chain.index:
                        # Only bond atoms in same chain 
                        newTopology.add_bond(prev_ca, new_ca)
                    prev_ca = new_ca
                atm_idx += 1

        self._ca_idxs = np.array(ca_idxs)
        self.topology = newTopology

    @property
    def top(self):
        return self.topology

    def map_traj(self, traj):
        """Create new trajectory"""
        ca_xyz = traj.xyz[:,self._ca_idxs[:,0],:]
        return md.Trajectory(ca_xyz, self.topology)

    @property
    def n_atomtypes(self):
        return len(self.atomtypes)

    def add_atoms(self, mass=1, radius=0.4, charge=0):
        name = "CA"
        mass = mass # amu
        radius = radius # nm
        charge = charge  # units??
        self.atoms = []
        self.atomtypes = []
        for atom in self.top.atoms:
            cg_atom = atom_types.CoarseGrainAtom(atom.index, name, 
                    atom.residue.index, atom.residue.name, radius, mass, charge)
            self.atoms.append(cg_atom)

            # Unique list of atom types.
            if cg_atom.name not in [ atm.name for atm in self.atomtypes ]:
                self.atomtypes.append(cg_atom)

    def _residue_to_atom_contacts(self, residue_contacts):
        atm_contacts = []
        for n,k in residue_contacts: 
            nres = self.topology.residue(n)
            kres = self.topology.residue(k)
            # Calpha-Calpha contact
            ca_n = nres.atom(0)
            ca_k = kres.atom(0)
            atm_contacts.append((ca_n, ca_k))
        return atm_contacts

    def _assign_sbm_angles(self):
        self._angles = []
        for chain in self.topology.chains:
            for i in range(chain.n_atoms - 2):
                self._angles.append(( chain.atom(i), chain.atom(i + 1), chain.atom(i + 2)))
        
    def _assign_sbm_dihedrals(self):
        self._improper_dihedrals = []
        self._dihedrals = []
        for chain in self.topology.chains:
            for i in range(chain.n_atoms - 3):
                self._dihedrals.append(( chain.atom(i), chain.atom(i + 1),\
                                   chain.atom(i + 2), chain.atom(i + 3)))

    def _assign_sbm_contacts(self, ref_traj_aa):
        residue_contacts = util.residue_contacts(ref_traj_aa)
        self._contact_pairs = self._residue_to_atom_contacts(residue_contacts)
    
    def _add_pairs(self, pairs):
        self._contact_pairs = []
        for n, k in pairs:
            natom = self.top.atom(n)
            katom = self.top.atom(k)
            self._contact_pairs.append([natom, katom])
    
    def add_disulfides(self, disulfides, simple=False):
        """ Add disulfide bonded interactions.
        
        Adds appropriate bond, angle and dihedral interactions.
        
        Args:
            disulfides (list): List of disulfide pairs between the 
                residues of the disulfides.
                
        """
        
        for pair in disulfides:
            res1 = self.top.residue(pair[0])
            res2 = self.top.residue(pair[1])
            respre = self.top.residue(pair[0]-1)
            respost = self.top.residue(pair[1]-1)
            
            #add c-alpha atoms
            cys1 = res1.atom(0)
            cys2 = res2.atom(0)
            #add c-beta
            ca1 = respre.atom(0)
            ca2 = respost.atom(0)
            
            #add bond between c-alpha of the Cys residues
            self.top.add_bond(cys1, cys2)
            
            if not simple:
                #add angular constraints between CYS;s and previous c-alphas
                self._angles.append((ca1, cys1, cys2))
                self._angles.append((cys1, cys2, ca2))
                
                #add dihedral constraints between CYS's and previous c-alphas
                self._dihedrals.append((ca1, cys1, cys2, ca2))
