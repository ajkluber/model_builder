import os
import numpy as np

import mdtraj as md

import viz_bonds

SUPPORTED_VERSIONS = ["4.5.4","4.5.4_sbm","4.6.5","4.6.5_sbm"]

class GromacsFiles(object):

    natively_supported_potentials = {"4.5.4":["LJ1210"]}

    def __init__(self, model, version=None):
        self.model = model
        self.version = version
        self.topfile = None

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
        #potentials that are supported when epsilon >=0
        self._switch_supported_pair_potentials = ["LJ12GAUSSIANTANH"]
        
    def _generate_index_file(self):
        top = self.model.mapping.topology
        self.index_ndx = "[ System ]\n"
        for i in range(top.n_atoms):
            atom = top.atom(i)
            self.index_ndx += "{:>4}".format(atom.index + 1)
            if ((i % 15) == 0) and (i > 0):
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
            if not self._check_supported(pot):
                self._tabled_pots.append(pot)

                table_name = "table_b{}.xvg".format(len(self._tabled_pots))
                self._tablenames.append(table_name)

                table = np.zeros((r.shape[0], 3), float)
                table[1:,0] = r[1:]
                table[10:,1] = pot.V(r[10:]) 
                table[10:,2] = -1.*pot.dVdr(r[10:]) 
                self._tables.append(table)
        self._n_tables = len(self._tablenames)
    
    def _check_supported(self, pot):
        """ Check if the potential requires a table file or not
        
        If the potential is supported, return True. No extra table file 
        will be written. If the potential is not supported, return 
        False. A table file of that potential will be generated in 
        _generate_interaction_tables(). Special cases may exist, such as 
        a gaussian-tanh function where it is supported for epsilon 
        greater than zero but not for less than zero. This is handled 
        inside nested if statements. Specific conditions can be added to 
        lists in the __init__ method.
        
        Parameters
        ----------
        pot : PairPotential
            Must have attribute `prefix_label` and `eps`.
        
        Returns
        -------
        supported : bool
            True if supported potential. False otherwise.
        
        """
        
        supported = True
        if pot.prefix_label in self._supported_pair_potentials:
            pass #unconditionally True            
        else:
            if pot.prefix_label in self._switch_supported_pair_potentials:
                if pot.eps >= 0:
                    pass #it is supported
                else:
                    supported = False
            else:
                supported = False
        
        return supported
        
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

    def write_simulation_files(self, path_to_tables=".", savetables=True, box_xyz=[20,20,20]):
        # Write the Hamiltonian Gromacs input file: topol.top
        if self.topfile is None:
            self.generate_topology()

        with open("topol.top", "w") as fout:
            fout.write(self.topfile)

        with open("index.ndx", "w") as fout:
            fout.write(self.index_ndx)

        if savetables:
            self._write_table_files(path_to_tables)

        # Save the starting configuration, but we have to fix the unitcell
        # information.
        if hasattr(self.model, "starting_traj"):
            self.model.starting_traj.save("conf.gro")
        else:
            self.model.ref_traj.save("conf.gro")

        with open("conf.gro", "r") as fin:
            temp = reduce(lambda x,y: x+y, fin.readlines()[:-1])
            temp += "{:>10f}{:>10f}{:>10f}\n".format(box_xyz[0], box_xyz[1], box_xyz[2])

        with open("conf.gro", "w") as fout:
            fout.write(temp)

        # files useful for visualizing
        self.model.ref_traj.save("ref.pdb")
        viz_bonds.write_bonds_tcl(self.model.mapping.top)

    def _write_table_files(self, path_to_tables):
        """Save table files"""

        cwd = os.getcwd()
        os.chdir(path_to_tables)
        if not os.path.exists("table.xvg"):
            np.savetxt("table.xvg", self.tablep, fmt="%16.15e")
        if not os.path.exists("tablep.xvg"):
            np.savetxt("tablep.xvg", self.tablep, fmt="%16.15e")

        for i in range(self._n_tables):
            if not os.path.exists(self._tablenames[i]):
                np.savetxt(self._tablenames[i], self._tables[i])
        os.chdir(cwd)

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
                    elif pot.prefix_label in ["LJ12GAUSSIAN", "LJ12GAUSSIANTANH"]:
                        func = 6
                        params = "{:>18.9e}{:>18.9e}{:>18.9e}{:>18.9e}".format(
                                    pot.eps, pot.r0, pot.width, pot.rNC**12)
                        assert pot.eps >= 0
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
