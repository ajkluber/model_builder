import numpy as np

import mdtraj as md
from mdtraj.core.topology import Topology

import util
import atom_types

class CalphaCbetaMapping(object):
    """Calpha Cbeta center-of-mass representation mapping"""

    def __init__(self, topology, use_chains=None):
        if use_chains is None:
            use_chains = range(len(topology._chains))

        self._ref_topology = topology.copy()

        # Build new topology
        newTopology = Topology()
        new_atm_idx = 0
        res_idx = 1
        prev_ca = None
        ca_idxs = []
        self._sidechain_idxs = []
        self._sidechain_mass = []
        self._chain_indices = []
        for chain_count,chain in enumerate(topology._chains):
            if chain_count in use_chains:
                newChain = newTopology.add_chain()
                for residue in chain._residues:
                    #resSeq = getattr(residue, 'resSeq', None) or residue.index
                    newResidue = newTopology.add_residue(residue.name, newChain, res_idx)
                    # map CA
                    new_ca = newTopology.add_atom('CA', md.core.element.get_by_symbol('C'),
                                        newResidue, serial=new_atm_idx)
                    self._chain_indices.append(chain_count)
                    if prev_ca is None:
                        prev_ca = new_ca
                    else:
                        # only bond atoms in the same chain.
                        if new_ca.residue.chain.index == prev_ca.residue.chain.index:
                            newTopology.add_bond(prev_ca, new_ca)
                        prev_ca = new_ca
                    try:
                        ca_idxs.append([[ atm.index for atm in residue.atoms if \
                                (atm.name == "CA") ][0], new_atm_idx ])
                    except:
                        print residue
                        print chain
                        for atm in residue.atoms:
                            atm.name
                        raise
                    new_atm_idx += 1

                    if residue.name == 'GLY':
                        self._sidechain_idxs.append([])
                        self._sidechain_mass.append([])
                    else:
                        # map CB
                        cb_name = "CB%s" % atom_types.residue_code[residue.name]
                        new_cb = newTopology.add_atom(cb_name, md.core.element.get_by_symbol('C'),
                                            newResidue, serial=new_atm_idx)
                        self._chain_indices.append(chain_count)

                        newTopology.add_bond(new_cb, new_ca)

                        self._sidechain_idxs.append([[ atm.index for atm in residue.atoms if \
                                    (atm.is_sidechain) and (atm.element.symbol != "H") ], new_atm_idx ])
                        self._sidechain_mass.append(np.array([ atm.element.mass for atm in residue.atoms if \
                                    (atm.is_sidechain) and (atm.element.symbol != "H") ]))
                        new_atm_idx += 1
                    res_idx += 1

        self._ca_idxs = np.array(ca_idxs)
        self.topology = newTopology
        assert self.topology.n_atoms == len(self._chain_indices)

    @property
    def top(self):
        return self.topology

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

    def _residue_to_atom_contacts(self, residue_contacts):
        atm_contacts = []
        for n,k in residue_contacts:
            nres = self.topology.residue(n)
            kres = self.topology.residue(k)
            # Calpha-Calpha contact
            ca_n = nres.atom(0)
            ca_k = kres.atom(0)
            atm_contacts.append([ca_n, ca_k])
            if (nres.name != "GLY") and (kres.name != "GLY"):
                # Cbeta-Cbeta contact
                cb_n = nres.atom(1)
                cb_k = kres.atom(1)
                atm_contacts.append([cb_n, cb_k])
        return atm_contacts

    def _assign_sbm_angles(self):

        self._angles = []
        chain_count = 0
        for chain in self.topology.chains:
            chain_count += 1
            # CA-CA-CA angles first
            ca_atoms = [ atom for atom in chain.atoms if atom.name == "CA" ]
            assert len(ca_atoms) == chain.n_residues
            for i in range(len(ca_atoms) - 2):
                self._angles.append((ca_atoms[i], ca_atoms[i + 1], ca_atoms[i + 2]))

            # CA-CA-CB and CB-CA-CA angles next
            for res_count,res in enumerate(chain.residues):
                if res.name != "GLY":
                    ca = res.atom(0)
                    cb = res.atom(1)
                    if res_count == 0:
                        # if terminal
                        next_ca = chain.residue(res_count + 1).atom(0)
                        self._angles.append((cb, ca, next_ca))
                    elif (res_count + 1) == chain.n_residues:
                        # if terminal
                        prev_ca = chain.residue(res_count - 1).atom(0)
                        self._angles.append((prev_ca, ca, cb))
                    else:
                        prev_ca = chain.residue(res_count - 1).atom(0)
                        next_ca = chain.residue(res_count + 1).atom(0)
                        self._angles.append((prev_ca, ca, cb))
                        self._angles.append((cb, ca, next_ca))

    def _assign_sbm_dihedrals(self):

        self._improper_dihedrals = []
        self._dihedrals = []
        for chain in self.topology.chains:
            #add the proper ca-ca-ca-ca dihedrals
            ca_atoms = [ atom for atom in chain.atoms if atom.name == "CA" ]
            for i in range(len(ca_atoms)-3):
                self._dihedrals.append((ca_atoms[i], ca_atoms[i+1], ca_atoms[i+2], ca_atoms[i+3]))
            #add improper dihedrals
            num_residues = chain.n_residues
            for res_count,res in enumerate(chain.residues):
                check = res_count == 0 #not first residue
                check = check or res_count+1 == num_residues #last residue
                check = check or res.name == "GLY" #GLY
                if not check:
                    idx = res_count
                    cj = chain.residue(idx-1).atom(0)
                    ck = chain.residue(idx+1).atom(0)
                    dih = (res.atom(0), cj, ck, res.atom(1))
                    self._improper_dihedrals.append(dih)

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

            #add c-alpha atoms
            ca1 = res1.atom(0)
            ca2 = res2.atom(0)
            #add c-beta
            cb1 = res1.atom(1)
            cb2 = res2.atom(1)

            #add bond between c-beta
            self.top.add_bond(cb1, cb2)

            #add angular constraints
            self._angles.append((ca1, cb1, cb2))
            self._angles.append((cb1, cb2, ca2))

            #add dihedral constraints
            self._dihedrals.append((ca1, cb1, cb2, ca2))

    @property
    def n_atomtypes(self):
        return len(self.atomtypes)

    def add_atoms(self):
        self.atoms = []
        self.atomtypes = []
        for chain in self.top._chains:
            for res in chain._residues:
                if res.name == "GLY":
                    cg_atom = atom_types.CoarseGrainAtom(res.atom(0).index, "CA",
                            res.index, res.name, 0.266, 1, 0)
                    self.atoms.append(cg_atom)
                else:
                    cg_atom = atom_types.CoarseGrainAtom(res.atom(0).index, "CA",
                            res.index, res.name, 0.266, 1, 0)
                    self.atoms.append(cg_atom)
                    radii = atom_types.residue_cacb_effective_interaction[res.name]
                    name = "CB%s" % atom_types.residue_code[res.name]
                    cg_atom = atom_types.CoarseGrainAtom(res.atom(1).index, name,
                            res.index, res.name, radii, 1, 0)
                    self.atoms.append(cg_atom)

        for cg_atom in self.atoms:
            # Unique list of atom types.
            if cg_atom.name not in [ atm.name for atm in self.atomtypes ]:
                self.atomtypes.append(cg_atom)

    def _atomidx_to_atom_contacts(self, pairs):
        atm_contacts = []
        for n,k in pairs:
            natom = self.top.atom(n)
            katom = self.top.atom(k)
            atm_contacts.append([natom, katom])
        return atm_contacts

    def _add_pairs(self, pairs):
        self._contact_pairs = self._atomidx_to_atom_contacts(pairs)

    def _add_atomtypes(self):
        pass
