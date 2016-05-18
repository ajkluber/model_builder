import numpy as np

import mdtraj as md
from mdtraj.core.topology import Topology

import util
import atom_types
        

class AwsemMapping(object):
    """Calpha Cbeta center-of-mass representation mapping"""

    def __init__(self, topology):
        self._ref_topology = topology.copy()

        # weights for creating glycine hydrogens (H-Beta, HB).
        # These weights are from Aram's scripts. 
        #self.HB_coeff_N = -0.946747
        #self.HB_coeff_CA = 2.50352
        #self.HB_coeff_C = -0.620388

        # There was something weird with the HB placement with Aram's
        # coefficients, so I have a simple solution of subtracting the center
        # of geometry of the N, CA, C atoms from the CA atom with some weight.
        self._HB_w = 3

        self._disulfides = []

        # Build new topology
        newTopology = Topology()
        CACBO_idxs = []
        HB_idxs = []
        res_charges = []
        res_idx = 1
        atm_idx = 0
        for chain in topology._chains:
            newChain = newTopology.add_chain()
            for residue in chain._residues:
                newTopology, atm_idx, res_idx = self._add_residue(newTopology,
                        newChain, residue, res_charges, res_idx, atm_idx,
                        CACBO_idxs, HB_idxs)

        self._charged_residues = res_charges
        self._HB_idxs = np.array(HB_idxs)
        self._CACBO_idxs = np.array(CACBO_idxs)
        self.topology = newTopology

    def _add_residue(self, newTopology, newChain, residue, res_charges,
                            res_idx, atm_idx, CACBO_idxs, HB_idxs):

        res_name = self._fix_nonstandard_resname(residue, res_charges)
        
        newResidue = newTopology.add_residue(res_name, newChain, res_idx)

        # Add atoms for residue 
        new_ca = newTopology.add_atom('CA', md.core.element.get_by_symbol('C'), newResidue, serial=atm_idx)
        CA_idx = [ atm.index for atm in residue.atoms if (atm.name == "CA") ][0]
        CACBO_idxs.append([CA_idx, new_ca.index])
        atm_idx += 1

        new_o = newTopology.add_atom('O', md.core.element.get_by_symbol('O'), newResidue, serial=atm_idx)
        O_idx = [ atm.index for atm in residue.atoms if (atm.name == "O") ][0]
        CACBO_idxs.append([O_idx, new_o.index])
        atm_idx += 1
        
        if residue.name == 'GLY':
            new_hb = newTopology.add_atom('HB', md.core.element.get_by_symbol('H'), newResidue, serial=atm_idx)
            N_idx = [ atm.index for atm in residue.atoms if (atm.name == "N") ][0]
            C_idx = [ atm.index for atm in residue.atoms if (atm.name == "C") ][0]
            HB_idxs.append([N_idx, CA_idx, C_idx, new_hb.index])
        else:
            new_cb = newTopology.add_atom('CB', md.core.element.get_by_symbol('C'), newResidue, serial=atm_idx)
            CB_idx = [ atm.index for atm in residue.atoms if (atm.name == "CB") ][0]
            CACBO_idxs.append([CB_idx, new_cb.index])

        # Add bonds to previous CA and O
        newTopology.add_bond(new_ca, new_o)

        if residue.name == "GLY":
            newTopology.add_bond(new_ca, new_hb)
        else:
            newTopology.add_bond(new_ca, new_cb)

        if residue.index > newChain.residue(0).index:
            prev_ca = [ atm for atm in newChain.residue(residue.index - 1 -\
                         newChain.residue(0).index).atoms if (atm.name == "CA") ][0]
            prev_o = [ atm for atm in newChain.residue(residue.index - 1 -\
                         newChain.residue(0).index).atoms if (atm.name == "O") ][0]
            newTopology.add_bond(prev_ca, new_ca)
            newTopology.add_bond(prev_o, new_ca)

        atm_idx += 1
        res_idx += 1

        return newTopology, atm_idx, res_idx

    def add_disulfides(self, disulfides):
        """Add disulfide bonds between corresponding CB atoms"""
        for pair in disulfides:
            res1 = self.top.residue(pair[0] - 1)
            res2 = self.top.residue(pair[1] - 1)

            cb1 = [ atom for atom in res1.atoms if (atom.name == "CB") ][0]
            cb2 = [ atom for atom in res2.atoms if (atom.name == "CB") ][0]

            self._disulfides.append([cb1, cb2])

    @property
    def top(self):
        return self.topology

    def map_traj(self, traj):
        """Return new Trajectory object with AWSEM topology and xyz"""
        # Direct slicing for CA, CB, O. 
        cacbo_xyz = np.zeros((traj.n_frames, self.topology.n_atoms, 3))
        cacbo_xyz[:, self._CACBO_idxs[:,1], :] = traj.xyz[:, self._CACBO_idxs[:,0], :]
        # HB is interpolated from other atoms.
        if not (len(self._HB_idxs) == 0):
            cacbo_xyz[:, self._HB_idxs[:,3], :] = (1. + (1. - (1./3))*self._HB_w)*traj.xyz[:, self._HB_idxs[:,1], :] -\
                                              (self._HB_w/3.)*traj.xyz[:, self._HB_idxs[:,0], :] -\
                                              (self._HB_w/3.)*traj.xyz[:, self._HB_idxs[:,2], :]

        return md.Trajectory(cacbo_xyz, self.top)

    def _fix_nonstandard_resname(self, residue, res_charges):
        """Convert non-standard residue names (that indicate charge state) to standard names"""
        res_name = residue.name
        new_res_name, charge = util._standard_resname_and_charge(res_name)
        if charge != 0:
            res_charges.append([residue.index + 1, charge])
        return new_res_name

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
                cb_n = nres.atom(2)
                cb_k = kres.atom(2)
                atm_contacts.append([cb_n, cb_k])
        return atm_contacts

    def mutate(self, mut_res_idx, mut_new_resname):
        """Mutate residue

        Parameters 
        ----------
        mut_res_idx : int
            Index of residue to mutate.
        mut_new_resname : str
            Three-letter code of residue to mutate to.
        """

        assert (self.topology.residue(mut_res_idx).name != mut_new_resname), "mutating the residue to itself!"

        # copy all residues except the mutated residue.
        # special cases: mutating to/from GLY.
        #   - mutating to GLY: need to rename CB -> HB 
        #   - mutating from GLY: rename HB -> CB. Also need to fix it's coordinates.

        # Build new topology
        newTopology = Topology()
        for chain in self.topology.chains:
            newChain = newTopology.add_chain()
            for residue in chain._residues:
                res_idx = residue.index
                if res_idx == mut_ref_idx:
                    # mutate residue
                    newResidue = newTopology.add_residue(mut_new_resname, newChain, res_idx)
                    for atom in residues.atoms:
                        if atom.name in ["CA", "O"]: 
                            # Copy over backbone atoms
                            newTopology.add_atom(atom.name,
                                    md.core.element.get_by_symbol(atom.element.symbol),
                                    newResidue, serial=atom.index)
                        else:
                            # Mutate sidechain atom if necessary
                            if (newResidue == "GLY") and (residue.name != "GLY"):
                                newTopology.add_atom("HB",
                                        md.core.element.get_by_symbol("H"),
                                        newResidue, serial=atom.index)
                            elif ((newResidue != "GLY") and (residue.name == "GLY")) or ((newResidue != "GLY") and (residue.name != "GLY")):
                                cb_atom = newTopology.add_atom("CB",
                                            md.core.element.get_by_symbol("C"),
                                            newResidue, serial=atom.index)
                                if ((newResidue != "GLY") and (residue.name == "GLY")):
                                    # TODO:
                                    # need new solution
                                    #HB_to_CB_idxs = [] 
                                    #N_idx = [ atm.index for atm in newResidue.atoms if (atm.name == "N") ][0]
                                    #C_idx = [ atm.index for atm in newResidue.atoms if (atm.name == "C") ][0]
                                    #HB_to_CB_idxs.append([N_ cb_atom.index])

                            else:
                                # mutating glycine to glycine. silly. 
                                newTopology.add_atom("HB",
                                        md.core.element.get_by_symbol("H"),
                                        newResidue, serial=atom.index)

                else:
                    # copy old residue atoms directly
                    newResidue = newTopology.add_residue(residue.name, newChain, res_idx)
                    for atom in residues.atoms:
                        newTopology.add_atom(atom.name, 
                                    md.core.element.get_by_symbol(atom.element.symbol), 
                                    newResidue, serial=atom.index)

        # The bond connectivity should stay identical
        for atm1, atm2 in self.topology._bonds:
            new_atm1 = newTopology.atom(atm1.index)
            new_atm2 = newTopology.atom(atm2.index)
            newTopology.add_bond(new_atm1, new_atm2)

        self._prev_topology = self.topology.copy()
        self.topology = newTopology

    def recenter_sidechains(self, traj, res_idxs):

        

class AwsemBackboneUnmapping(object):
    """Calpha Cbeta center-of-mass representation mapping"""

    def __init__(self, topology):
        self._ref_topology = topology.copy()

        self._N_coeff = np.array([0.48318, 0.70328, -0.18643])  # CA_i-1, CA_i, O_i-1
        self._C_coeff = np.array([0.44365, 0.23520, 0.32115])   # CA_i, CA_i+1, O_i
        self._H_coeff = np.array([0.84100, 0.89296, -0.73389])  # CA_i-1, CA_i, O_i-1

        newTopology = Topology()
        CACBO_idxs = []; N_idxs = []; C_idxs = []; H_idxs = []
        res_idx = 1
        atm_idx = 0
        prev_ca = None
        prev_o = None
        for chain in topology._chains:
            newChain = newTopology.add_chain()
            for residue in chain._residues:
                newTopology, atm_idx, res_idx = self._add_residue(newTopology, 
                        newChain, residue, chain, res_idx, atm_idx, 
                        N_idxs, C_idxs, H_idxs, CACBO_idxs, prev_ca, prev_o)

        self._CACBO_idxs = np.array(CACBO_idxs)
        self._N_idxs = np.array(N_idxs)
        self._C_idxs = np.array(C_idxs)
        self._H_idxs = np.array(H_idxs)
        self.topology = newTopology

    def _add_residue(self, newTopology, newChain, residue, chain, 
                res_idx, atm_idx, N_idxs, C_idxs, H_idxs, CACBO_idxs, prev_ca, prev_o):

        res_name = self._fix_nonstandard_resname(residue, [])
        
        newResidue = newTopology.add_residue(res_name, newChain, res_idx)

        # Add atoms that are directly mapped 
        new_ca = newTopology.add_atom('CA', md.core.element.get_by_symbol('C'), 
                                    newResidue, serial=atm_idx)
        CA_idx = [ atm.index for atm in residue.atoms if (atm.name == "CA") ][0]
        CACBO_idxs.append([CA_idx, new_ca.index])
        atm_idx += 1

        new_o = newTopology.add_atom('O', md.core.element.get_by_symbol('O'), 
                                    newResidue, serial=atm_idx)
        O_idx = [ atm.index for atm in residue.atoms if (atm.name == "O") ][0]
        CACBO_idxs.append([O_idx, new_o.index])
        newTopology.add_bond(new_ca, new_o)
        atm_idx += 1

        if residue.index > chain.residue(0).index:
            # add N atom if not the N-terminus 
            prev_res = chain.residue(residue.index - 1)
            new_n = newTopology.add_atom('N', md.core.element.get_by_symbol('N'), 
                                        newResidue, serial=atm_idx)
            
            prev_CA_idx = [ atm.index for atm in prev_res.atoms if (atm.name == "CA") ][0]
            prev_O_idx = [ atm.index for atm in prev_res.atoms if (atm.name == "O") ][0]
            N_idxs.append([prev_CA_idx, CA_idx, prev_O_idx, new_n.index])

            prev_c = [ atm for atm in newTopology.residue(newResidue.index - 1).atoms if (atm.name == "C") ][0]
            newTopology.add_bond(new_n, new_ca)
            newTopology.add_bond(new_n, prev_c)
            atm_idx += 1

            new_h = newTopology.add_atom('H', md.core.element.get_by_symbol('H'), 
                                        newResidue, serial=atm_idx)
            H_idxs.append([prev_CA_idx, CA_idx, prev_O_idx, new_h.index])
            newTopology.add_bond(new_n, new_h)

            atm_idx += 1

        if residue.index < chain.residue(chain.n_residues - 1).index:
            # add C atom if not the C-terminus
            next_res = chain.residue(residue.index + 1)
            new_c = newTopology.add_atom('C', md.core.element.get_by_symbol('C'), 
                                        newResidue, serial=atm_idx)
            next_CA_idx = [ atm.index for atm in next_res.atoms if (atm.name == "CB") ][0]
            C_idxs.append([CA_idx, next_CA_idx, O_idx, new_c.index])
            newTopology.add_bond(new_c, new_ca)
            newTopology.add_bond(new_c, new_o)
            atm_idx += 1

        if residue.name == 'GLY':
            new_hb = newTopology.add_atom('HB', md.core.element.get_by_symbol('H'), 
                                        newResidue, serial=atm_idx)
            H_idx = [ atm.index for atm in residue.atoms if (atm.name == "HB") ][0]
            CACBO_idxs.append([H_idx, new_hb.index])
            newTopology.add_bond(new_ca, new_hb)
            atm_idx += 1
        else:
            new_cb = newTopology.add_atom('CB', md.core.element.get_by_symbol('C'), 
                                        newResidue, serial=atm_idx)
            CB_idx = [ atm.index for atm in residue.atoms if (atm.name == "CB") ][0]
            CACBO_idxs.append([CB_idx, new_cb.index])
            newTopology.add_bond(new_ca, new_cb)
            atm_idx += 1

        # Add bonds to previous CA and O
        if (prev_ca is None) and (prev_o is None):
            prev_ca = new_ca
            prev_o = new_o
        else:
            # Only bond atoms in the same chain
            if prev_ca.residue.chain.index == new_ca.residue.chain.index:
                newTopology.add_bond(prev_ca, new_ca)
            if prev_o.residue.chain.index == new_ca.residue.chain.index:
                newTopology.add_bond(prev_o, new_ca)
            prev_ca = new_ca
            prev_o = new_o

        atm_idx += 1
        res_idx += 1

        return newTopology, atm_idx, res_idx

    def _fix_nonstandard_resname(self, residue, res_charges):
        """Convert non-standard residue names (that indicate charge state) to standard names"""
        res_name = residue.name
        new_res_name, charge = util._standard_resname_and_charge(res_name)
        if charge != 0:
            res_charges.append([residue.index + 1, charge])
        return new_res_name

    def map_traj(self, traj):
        """Map coarse-grained coordinates to all-atom backbone coordinates"""
        newxyz = np.zeros((traj.n_frames, self.topology.n_atoms, 3))
        newxyz[:,self._CACBO_idxs[:,1],:] = traj.xyz[:,self._CACBO_idxs[:,0]]

        # interpolate N, C, and H
        newxyz[:, self._N_idxs[:,3], :] = self._N_coeff[0]*traj.xyz[:, self._N_idxs[:,0], :] +\
                                        self._N_coeff[1]*traj.xyz[:, self._N_idxs[:,1], :] +\
                                        self._N_coeff[2]*traj.xyz[:, self._N_idxs[:,2], :]
        newxyz[:, self._C_idxs[:,3], :] = self._C_coeff[0]*traj.xyz[:, self._C_idxs[:,0], :] +\
                                        self._C_coeff[1]*traj.xyz[:, self._C_idxs[:,1], :] +\
                                        self._C_coeff[2]*traj.xyz[:, self._C_idxs[:,2], :]
        newxyz[:, self._H_idxs[:,3], :] = self._H_coeff[0]*traj.xyz[:, self._H_idxs[:,0], :] +\
                                        self._H_coeff[1]*traj.xyz[:, self._H_idxs[:,1], :] +\
                                        self._H_coeff[2]*traj.xyz[:, self._H_idxs[:,2], :]
        return md.Trajectory(newxyz, self.topology)

