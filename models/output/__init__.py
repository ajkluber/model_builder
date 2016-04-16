import os
import numpy as np

SUPPORTED_VERSIONS = ["4.5.4","4.5.4_sbm","4.6.5","4.6.5_sbm"]

class GromacsFiles(object):

    natively_supported_potentials = {"4.5.4":["LJ1210",]}

    def __init__(self, model, version=None):
        self.model = model
        self.version = version

        # check compatibility of interactions with this version
        # of gromacs

        # check for interactions that need to be tabled.

        # determine the right lookup numbers for each interactions 
        self._bond_funcs = {"HARMONIC_BOND":1}
        self._angle_funcs = {"HARMONIC_ANGLE":1}
        self._dihedral_funcs = {"COSINE_DIHEDRAL":1,
                                "HARMONIC_DIHEDRAL":2}

        self._supported_pair_potentials = ["LJ12", "LJ1210",
                                            "GAUSSIAN", "LJ12GAUSSIAN"]

    def _generate_index_file(self):
        top = self.model.mapping.topology
        self.index_ndx = "[ System ]\n"
        for i in range(top.n_atoms):
            atom = top.atom(i)
            self.index_ndx += "{:>4}".format(atom.index + 1)
            if (i % 15) == 0:
                self.index_ndx += "\n" 
        self.index_ndx += "\n" 

    def _generate_interaction_tables(self):
        """Generates tables of user-defined potentials"""

        # Determine which interactions need to be tabled 
        self.tablep = self._get_LJ1210_table()
        self.LJtable = self.tablep

        r = np.arange(0, 20.0, 0.002)
        self._tabled_pots = []
        self._tables = []
        self._tablenames = []
        for i in range(self.model.Hamiltonian.n_pairs):
            pot = self.model.Hamiltonian._pairs[i]
            if pot.prefix_label not in self._supported_pair_potentials:
                self._tabled_pots.append(pot)

                table_name = "table_b{}.xvg".format(len(self._tabled_pots))
                self._tablenames.append(table_name)

                table = np.zeros((r.shape[0], 3), float)
                table[1:,0] = r[1:]
                table[10:,1] = pot.V(r[10:]) 
                table[10:,2] = -1.*pot.dVdr(r[10:]) 
                self._tables.append(table)
        self._n_tables = len(self._tablenames)

    def _get_LJ1210_table(self):
        """ LJ1210 interaction table """ 
        r = np.arange(0.0, 20.0, 0.002)
        r[0] = 1
        table = np.zeros((len(r),7),float)
        table[:,0] = r
        table[:,1] = 1./r
        table[:,2] = 1./(r**2)
        table[:,3] = -1./(r**10)
        table[:,4] = -10./(r**11)
        table[:,5] = 1./(r**12)
        table[:,6] = 12./(r**13)
        table[:5,1:] = 0
        table[0,0] = 0
        return table

    def write_simulation_files(self, path_to_tables="."):
        # Write the Hamiltonian Gromacs input file: topol.top
        self.generate_topology()

        with open("topol.top", "w") as fout:
            fout.write(self.topfile)

        with open("index.ndx", "w") as fout:
            fout.write(self.index_ndx)

        self._write_table_files(path_to_tables)

        # Save the starting configuration, but we have to fix the unitcell
        # information.
        self.model.ref_traj.save("conf.gro")
        with open("conf.gro", "r") as fin:
            temp = reduce(lambda x,y: x+y, fin.readlines()[:-1])
            temp += "{:>10f}{:>10f}{:>10f}\n".format(20,20,20)

        with open("conf.gro", "w") as fout:
            fout.write(temp)

    def _write_table_files(self, path_to_tables):
        """Save table files"""

        if not os.path.exists("{}/table.xvg".format(path_to_tables)):
            np.savetxt("{}/table.xvg".format(path_to_tables), self.tablep, fmt="%16.15e")
        if not os.path.exists("{}/tablep.xvg".format(path_to_tables)):
            np.savetxt("{}/tablep.xvg".format(path_to_tables), self.tablep, fmt="%16.15e")

        for i in range(self._n_tables):
            if not os.path.exists("{}/{}".format(path_to_tables,self._tablenames[i])):
                np.savetxt("{}/{}".format(path_to_tables,self._tablenames[i]), self._tables[i])

    def _get_atomtypes_top(self):
        """ Generate the [ atoms ] top."""
        atomtypes_top = " [ atomtypes ]\n"
        atomtypes_top += " ;name   mass  charge ptype      c6               c12\n"
        for atomtype in self.model.mapping.atomtypes:
            atomtypes_top += "{:<3}  {:>8.3f}{:>8.3f} {} {:>18.9e}{:>18.9e}\n".format(
                            atomtype.name, atomtype.mass, atomtype.charge,
                            atomtype.ptype, atomtype.c6, atomtype.c12)
        return atomtypes_top

    def _get_atoms_top(self):
        """ Generate the [ atoms ] top."""
        atoms_top = " [ atoms ]\n"
        atoms_top += " ;  nr  type resnr  res atom   cgnr  charge    mass\n"
        for atom in self.model.mapping.atoms:
            atoms_top += " {:>5d}{:>4}{:>8d}{:>5}{:>4}{:>8}{:>8.3f}{:>8.3f}\n".format(
                                atom.index + 1, atom.name, atom.residx + 1, 
                                atom.resname, atom.name, atom.index + 1, atom.charge, atom.mass)
        return atoms_top

    def _get_bonds_top(self):
        """ Generate the [ bonds ] top."""
        if self.model.Hamiltonian.n_bonds == 0:
            bonds_top = ""
        else:
            bonds_top = " [ bonds ]\n"
            bonds_top += " ;  ai     aj func     r0(nm)            Kb\n"
            for bond in self.model.Hamiltonian.bonds:
                func = self._bond_funcs[bond.prefix_label]
                bonds_top += "{:>6} {:>6}{:>2}{:>18.9e}{:>18.9e}\n".format(
                                    bond.atmi.index + 1, bond.atmj.index + 1,
                                    func, bond.r0, bond.kb) 
        return bonds_top

    def _get_tabled_top(self):
        """ Generate the topology files to specify table interactions. """
        # Add special nonbonded table interactions. 
        if self._n_tables == 0:
            tabled_string = ""
        else:
            tabled_string = "; tabled interactions pairs below\n"
            for i in range(self._n_tables):
                pot = self._tabled_pots[i]
                tabled_string += "{:>6} {:>6}{:>2}{:>18}{:>18.9e}\n".format(
                              pot.atmi.index + 1, pot.atmj.index + 1, 9, i + 1, 1)
        return tabled_string

    def _get_angles_top(self):
        """ Generate the [ angles ] top."""
        if self.model.Hamiltonian.n_angles == 0:
            angles_top = ""
        else:
            angles_top = " [ angles ]\n"
            angles_top += " ;  ai     aj     ak func     th0(deg)            Ka\n"
            for angle in self.model.Hamiltonian.angles:
                if angle.prefix_label == "HARMONIC_ANGLE":
                    ka = angle.ka
                    func = self._angle_funcs[angle.prefix_label]
                    angles_top += "{:>6} {:>6} {:>6}{:>2}{:>18.9e}{:>18.9e}\n".format(
                                    angle.atmi.index + 1, angle.atmj.index + 1, 
                                    angle.atmk.index + 1, func, 
                                    angle.theta0*(180./np.pi), ka)
        return angles_top

    def _get_dihedrals_top(self):
        """ Generate the [ dihedrals ] top."""
        if self.model.Hamiltonian.n_dihedrals == 0:
            dihedrals_top = ""
        else:
            dihedrals_top = " [ dihedrals ]\n"
            dihedrals_top += " ;  ai     aj     ak     al func    phi0(deg)           Kd        (mult)\n"
            for dih in self.model.Hamiltonian.dihedrals:
                func = self._dihedral_funcs[dih.prefix_label]
                if dih.prefix_label == "COSINE_DIHEDRAL":
                    #first convert to degrees, then Gromacs angles, then multiplicity
                    phi_s = dih.mult*(180. + dih.phi0*(180./np.pi)) #first convert to Gromacs angles, then multiply by multiplicity
                    
                    dihedrals_top += "{:>6} {:>6} {:>6} {:>6}{:>2}{:>18.9e}{:>18.9e}{:>3d}\n".format(
                                    dih.atmi.index + 1, dih.atmj.index + 1,
                                    dih.atmk.index + 1, dih.atml.index + 1,
                                    func, phi_s, dih.kd, dih.mult)
                elif dih.prefix_label == "HARMONIC_DIHEDRAL":
                    phi_s = dih.phi0*(180./np.pi)
                    #kd = dih.kd*((np.pi/180.)**2) 
                    kd = dih.kd
                    dihedrals_top += "{:>6} {:>6} {:>6} {:>6}{:>2}{:>18.9e}{:>18.9e}\n".format(
                                    dih.atmi.index + 1, dih.atmj.index + 1,
                                    dih.atmk.index + 1, dih.atml.index + 1,
                                    func, phi_s, kd)
                else:
                    print "Warning: unknown dihedral interaction for: {}".format(dih.describe())
        return dihedrals_top

    def _get_pairs_top(self):
        """ Get the [ pairs ] top for SBM Gromacs """
        if self.model.Hamiltonian.n_pairs == 0:
            pairs_top = ""
        else:
            pairs_top = " [ pairs ]\n"
            pairs_top += " ;   i      j type      c10               c12  \n"
            for i in range(self.model.Hamiltonian.n_pairs):
                pot = self.model.Hamiltonian._pairs[i]
                # How are LJ126 interactions defined? By atomtypes(?)
                if pot not in self._tabled_pots:
                    atm_idxs = "{:>6} {:>6}".format(pot.atmi.index + 1, pot.atmj.index + 1)
                    if pot.prefix_label == "LJ1210":
                        func = 1
                        c12 = pot.eps*5.0*(pot.r0**12)
                        c10 = pot.eps*6.0*(pot.r0**10)
                        params = "{:>18.9e}{:>18.9e}".format(c10, c12)
                    elif pot.prefix_label == "LJ12":
                        func = 1
                        c12 = pot.eps*(pot.r0**12)
                        c10 = 0
                        params = "{:>18.9e}{:>18.9e}".format(c10, c12)
                    elif pot.prefix_label == "GAUSSIAN":
                        func = 5
                        params = "{:>18.9e}{:>18.9e}{:>18.9e}".format(
                                    pot.eps, pot.r0, pot.width)
                    elif pot.prefix_label == "LJ12GAUSSIAN":
                        func = 6
                        params = "{:>18.9e}{:>18.9e}{:>18.9e}{:>18.9e}".format(
                                    pot.eps, pot.r0, pot.width, pot.rNC**12)
                    else:
                        print "Warning: interaction is not supported: {}".format(pot.describe())
                    pairs_top += "{}{:>2}{}\n".format(atm_idxs, func, params)
        return pairs_top
    
    def _get_exclusions_top(self):
        """ Get [ exclusions ] top""" 
        if self.model.Hamiltonian.n_pairs == 0:
            exclusions_top = ""
        else:
            exclusions_top = " [ exclusions ]\n"
            for pot in self.model.Hamiltonian.pairs:
                exclusions_top += "{:>6} {:>6}\n".format(pot.atmi.index + 1, 
                                                         pot.atmj.index + 1)
        return exclusions_top

    def generate_topology(self):
        """ Return a structure-based topology file. Currently only for one molecule. """

        self._generate_interaction_tables()
        self._generate_index_file()

        top =  " ; Structure-based  topology file for Gromacs:\n"
        top += " [ defaults ]\n"
        top += " ;nbfunc    comb-rule gen-pairs\n"
        top += "      1           1 no\n\n"
        top += "{}\n".format(self._get_atomtypes_top())
        top += " [ moleculetype ]\n"
        top += " ;name                 nrexcl\n"
        top += " Macromolecule           3\n\n"

        top += "{}\n".format(self._get_atoms_top())
        top += "{}".format(self._get_bonds_top())
        top += "{}\n".format(self._get_tabled_top())
        top += "{}\n".format(self._get_angles_top())
        top += "{}\n".format(self._get_dihedrals_top())
        top += "{}\n".format(self._get_pairs_top())
        top += "{}\n".format(self._get_exclusions_top())

        top += " [ system ]\n"
        top += " ;   name\n"
        top += " Macromolecule\n\n"
        top += " [ molecules ]\n"
        top += " ;   name      n_molec \n"
        top += " Macromolecule 1\n\n" # Need to account for multiple molecules (?)

        self.topfile = top


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







