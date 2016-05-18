import os
import numpy as np

import mdtraj as md

class AWSEMLammpsFiles(object):

    def __init__(self, model, version=""):
        self.model = model
        self.version = version

        self._bond_coeffs = [["harmonic", 200.,3.77], 
                            ["harmonic", 200.,2.41],
                            ["harmonic", 200.,2.50],
                            ["harmonic", 200.,1.54],
                            ["harmonic", 200.,1.54]]

        self._disulfide_bond_coeffs = ["harmonic", 200., 3.72]

        self._bond_type_idxs = {"CACA":1, "CAO":2, "OCA":3, 
                                "CACB":4, "CBCA":4, "CAHB":5, "HBCA":5,
                                "disulfide":6}


        self._masses = [27., 14., 28., 60., 2.]

        self._atom_type_idxs = {"CA":1, "N":2, "O":3, "CB":4, "HB":5}

    def _get_header_string(self, box_xyz_lo_hi):
        top = self.model.mapping.top

        n_disulfides = len(self.model.mapping._disulfides)

        header_string = "LAMMPS written by model_builder\n\n"
        header_string += "{:>12d}  atoms\n".format(top.n_atoms)
        header_string += "{:>12d}  bonds\n".format(top.n_bonds + n_disulfides)
        header_string += "{:>12d}  angles\n".format(0)
        header_string += "{:>12d}  dihedrals\n".format(0)
        header_string += "{:>12d}  impropers\n\n".format(0)

        header_string += "{:>12d}  atom types\n".format(5)
        if n_disulfides > 0:
            header_string += "{:>12d}  bond types\n".format(6)
        else:
            header_string += "{:>12d}  bond types\n".format(5)
        header_string += "{:>12d}  angle types\n".format(0)
        header_string += "{:>12d}  dihedral types\n".format(0)
        header_string += "{:>12d}  improper types\n\n".format(0)

        header_string += "{:<8.1f} {:<8.1f}xlo xhi\n".format(box_xyz_lo_hi[0][0],box_xyz_lo_hi[0][1])
        header_string += "{:<8.1f} {:<8.1f}ylo yhi\n".format(box_xyz_lo_hi[1][0],box_xyz_lo_hi[1][1])
        header_string += "{:<8.1f} {:<8.1f}zlo zhi\n\n".format(box_xyz_lo_hi[2][0],box_xyz_lo_hi[2][1])
        return header_string

    def _get_masses_string(self):
        masses_string = "Masses\n\n"
        for i in range(len(self._masses)):
            masses_string += "{:>12d}    {:.2f}\n".format(i + 1, self._masses[i])
        masses_string += "\n"
        return masses_string

    def _get_atoms_string(self):
        if hasattr(self.model, "starting_traj"):
            xyz = 10.*self.model.starting_traj.xyz[0]
        elif hasattr(self.model, "ref_traj"):
            xyz = 10.*self.model.ref_traj.xyz[0]
        else:
            raise AttributeError("need to set intial conditions (ref_traj or starting_traj) to write")

        charge = 0
        atoms_string = "Atoms\n\n"
        for i in range(self.model.mapping.top.n_atoms):
            atom = self.model.mapping.top.atom(i)
            atom_type = self._atom_type_idxs[atom.name]
            atoms_string += "{:>12d}{:>5d}{:>5d}{:>5d} {:>6f} {:>10f} {:>10f} {:>10f}\n".format(
                        i + 1, atom.residue.chain.index + 1, atom.residue.index + 1, 
                        atom_type, charge, xyz[i, 0], xyz[i, 1], xyz[i, 2])
        atoms_string += "\n"
        return atoms_string

    def _get_bond_coeffs_string(self):
        bond_coeff_string = "Bond Coeffs\n\n"
        for i in range(len(self._bond_coeffs)):
            bond_coeff_string += "{:>12d} {:>6.3f} {:>6.3f}\n".format(i + 1, 
                   self._bond_coeffs[i][1], self._bond_coeffs[i][2])

        if len(self.model.mapping._disulfides) > 0:
            bond_coeff_string += "{:>12d} {:>6.3f} {:>6.3f}\n".format(i + 2, 
                    self._disulfide_bond_coeffs[1], self._disulfide_bond_coeffs[2])

        bond_coeff_string += "\n"
        return bond_coeff_string

    def _get_bonds_string(self):
        bond_string = "Bonds\n\n"
        for i in range(len(self.model.mapping.top._bonds)):
            atmi, atmj = self.model.mapping.top._bonds[i]
            bond_type = self._bond_type_idxs[atmi.name + atmj.name]
            bond_string += "{:>11d} {:>5d} {:>5d} {:>5d}\n".format(i + 1, 
                            bond_type, atmi.index + 1, atmj.index + 1)

        if len(self.model.mapping._disulfides) > 0:
            start_idx = i + 2
            for i in range(len(self.model.mapping._disulfides)):
                atmi, atmj = self.model.mapping._disulfides[i]
                bond_type = self._bond_type_idxs["disulfide"]
                bond_string += "{:>11d} {:>5d} {:>5d} {:>5d}\n".format(
                        start_idx + i, bond_type, atmi.index + 1, atmj.index + 1)

        bond_string += "\n"
        return bond_string

    def generate_topology(self, box_xyz_lo_hi=[[-200,200],[-200,200],[-200,200]]):
        """Generate topology file for AWSEM-LAMMPS simulation"""
        
        top_string = self._get_header_string(box_xyz_lo_hi)
        top_string += self._get_masses_string()
        top_string += self._get_atoms_string()
        top_string += self._get_bond_coeffs_string()
        top_string += self._get_bonds_string()

        self.topfile = top_string

    def write_simulation_files(self, ref_traj_aa, topfilename, seqfilename):

        self.generate_topology()

        if hasattr(self.model, "ref_traj"):
            traj = self.model.ref_traj
        elif hasattr(self.model, "starting_traj"):
            traj = self.model.starting_traj
        else:
            raise AttributeError("need to set intial conditions (ref_traj or starting_traj) to write")

        fasta = traj.top.to_fasta()

        with open("{}".format(seqfilename),"w") as fout:
            for line in fasta:
                fout.write("{}\n".format(line))

        with open("charge_on_residues.dat", "w") as fout:
            fout.write("{:d}\n".format(len(self.model.mapping._charged_residues)))
            for res in self.model.mapping._charged_residues:
                fout.write("{:6d}   {:8.4f}\n".format(res[0], res[1]))

        # compute secondary structure from a reference structure
        dssp = ("".join(md.compute_dssp(ref_traj_aa)[0])).replace("C","-")
        assert len(dssp) == sum([ len(x) for x in fasta ]), "Number of residues in reference different than expected"
        with open("ssweight", "w") as fout:
            for ss in dssp:
                if ss == "H": 
                    helix = 1.
                    sheet = 0.
                elif ss == "E":
                    helix = 0.
                    sheet = 1.
                else:
                    helix = 0.
                    sheet = 0.
                fout.write("{:.1f} {:.1f}\n".format(helix, sheet))

        with open("jpred", "w") as fout:
            start = 0
            for i in range(len(fasta)):
                chain_length = len(fasta[i]) 
                fout.write("{}\n".format(fasta[i]))
                fout.write("{}\n".format(dssp[start:start+chain_length]))
                start += chain_length


        with open("{}".format(topfilename),"w") as fout:
            fout.write(self.topfile)
