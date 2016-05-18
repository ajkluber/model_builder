import os
import numpy as np

import mdtraj as md

class LammpsFiles(object):

    def __init__(self, model, version=""):
        self.model = model
        self.version = version

        # determine the right lookup numbers for each interactions 
        #self._bond_funcs = {"HARMONIC_BOND":1}
        #self._angle_funcs = {"HARMONIC_ANGLE":1}
        #self._dihedral_funcs = {"COSINE_DIHEDRAL":1,
        #                        "HARMONIC_DIHEDRAL":2}

        #self._supported_pair_potentials = ["LJ12", "LJ1210",
        #                                    "GAUSSIAN", "LJ12GAUSSIAN"]

    def generate_topology(self):

        mapping = self.model.mapping
        top = self.model.mapping.top
        if hasattr(self.model, "starting_traj"):
            xyz = self.model.starting_traj.xyz[0]
        elif hasattr(self.model, "ref_traj"):
            xyz = self.model.ref_traj.xyz[0]
        else:
            raise AttributeError("need to set intial conditions (ref_traj or starting_traj) to write")
        
        top_string = "LAMMPS written by model_builder\n\n"
        top_string += "{:>12d}  atoms\n".format(top.n_atoms)
        top_string += "{:>12d}  bonds\n".format(top.n_bonds)
        top_string += "{:>12d}  angles\n".format(0)
        top_string += "{:>12d}  dihedrals\n".format(0)
        top_string += "{:>12d}  impropers\n\n".format(0)

        top_string += "{:>12d}  atom types\n".format(mapping.n_atomtypes)
        top_string += "{:>12d}  bond types\n\n".format(top.n_bonds)

        top_string += "{:<8.1f} {:<8.1f}xlo xhi\n".format(-200., 200.)
        top_string += "{:<8.1f} {:<8.1f}ylo yhi\n".format(-200., 200.)
        top_string += "{:<8.1f} {:<8.1f}zlo zhi\n\n".format(-200., 200.)

        top_string += "Masses\n\n"
        atomtype_map = {}
        for i in range(mapping.n_atomtypes):
            atomtype_map[mapping.atomtypes[i].name] = i + 1
            top_string += "{:>12d}    {:.2f}\n".format(i + 1, mapping.atomtypes[i].mass)
        
        top_string += "\n"
        top_string += "Atoms\n\n"
        for i in range(top.n_atoms):
            atom = top.atom(i)
            atom_type = mapping.atoms[i]
            mass = atom_type.mass
            charge = atom_type.charge
            type_idx = atomtype_map[atom_type.name]
            top_string += "{:>12d}{:>5d}{:>5d}{:>5d} {:>6f} {:>10f} {:>10f} {:>10f}\n".format(
                        i + 1, atom.residue.chain.index + 1, atom.residue.index + 1, 
                        type_idx, charge, xyz[i,0], xyz[i,1], xyz[i,2])
        top_string += "\n"

        top_string += "Bond Coeffs\n\n"
        for i in range(self.model.Hamiltonian.n_bonds):
            # CHECK UNITS
            bond = self.model.Hamiltonian._bonds[i]
            top_string += "{:>12d} {:>6.3f} {:>6.3f}\n".format(i + 1, bond.kb, bond.r0*10.)
        top_string += "\n"

        top_string += "Bonds\n\n"
        for i in range(self.model.Hamiltonian.n_bonds):
            bond = self.model.Hamiltonian._bonds[i]
            top_string += "{:>11d} {:>5d} {:>5d} {:>5d}\n".format(i + 1, i + 1, bond.atmi.index + 1, bond.atmj.index + 1)
        top_string += "\n"

        self.topfile = top_string



