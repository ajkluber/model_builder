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

                table_name = "table_b{}.xvg".format(len(self._tabled_pots) + 1)
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

    def write_simulation_files(self):
        pass

    def _get_atomtypes_top(self):
        """ Generate the [ atoms ] top."""
        atomtypes_top = " [ atomtypes ]\n"
        atomtypes_top += " ;name   mass  charge ptype      c6               c12\n"
        for atomtype in self.model.mapping.atomtypes:
            atomtypes_top += " {}  {:>8.3f}{:>8.3f} {} {:>18.9e}{:>18.9e}\n".format(
                            atomtype.name, atomtype.mass, atomtype.charge,
                            atomtype.ptype, atomtype.c6, atomtype.c12)
        return atomtypes_top

    def _get_atoms_top(self):
        """ Generate the [ atoms ] top."""
        top = self.model.mapping.topology
        atoms_top = " [ atoms ]\n"
        atoms_top += " ;  nr  type resnr  res atom   cgnr  charge    mass\n"
        for atom in top.atoms:
            atoms_top += " {:>5d}{:>4}{:>8d}{:>5}{:>4}{:>8}{:>8.3f}{:>8.3f}\n".format(
                                atom.index + 1, atom.name, atom.residue.index + 1, 
                                atom.residue.name, atom.name, atom.index + 1, 0.0, 1.0)
        return atoms_top

    def _get_bonds_top(self):
        """ Generate the [ bonds ] top."""
        bonds_top = " [ bonds ]\n"
        bonds_top += " ;  ai     aj func     r0(nm)            Kb\n"
        for bond in self.model.Hamiltonian.bonds:
            func = self._bond_funcs[bond.prefix_label]
            bonds_top += "{:>6} {:>6}{:>2}{:>18.9e}{:>18.9e}\n".format(
                                bond.atmi.index + 1, bond.atmj.index + 1, func, bond.r0, bond.kb) 
        return bonds_top

    def _get_tabled_top(self):
        """ Generate the topology files to specify table interactions. """
        # Add special nonbonded table interactions. 
        if self._n_tables != 0:
            tabled_string = "; tabled interactions pairs below\n"
            for i in range(self._n_tables):
                pot = self._tabled_pots[i]
                tabled_string += "{:>6} {:>6}{:>2}{:>18}{:>18.9e}\n".format(
                              pot.atmi.index + 1, pot.atmj.index + 1, 9, i + 1, 1)
        else:
            tabled_string = ""
        return tabled_string

    def _get_angles_top(self):
        """ Generate the [ angles ] top."""
        angles_top = " [ angles ]\n"
        angles_top += " ;  ai     aj     ak func     th0(deg)            Ka\n"
        for angle in self.model.Hamiltonian.angles:
            func = self._angle_funcs[angle.prefix_label]
            angles_top += "{:>6} {:>6} {:>6}{:>2}{:>18.9e}{:>18.9e}\n".format(
                            angle.atmi.index + 1, angle.atmj.index + 1, 
                            angle.atmk.index + 1, func, 
                            angle.theta0, angle.ka)
        return angles_top

    def _get_dihedrals_top(self):
        """ Generate the [ dihedrals ] top."""
        dihedrals_top = " [ dihedrals ]\n"
        dihedrals_top += " ;  ai     aj     ak     al func    phi0(deg)           Kd        (mult)\n"
        for dih in self.model.Hamiltonian.dihedrals:
            func = self._dihedral_funcs[dih.prefix_label]
            #self._dihedrals_top += "%6d %6d %6d %6d%2d%18.9e%18.9e%2d\n" %  \
            if dih.prefix_label == "COSINE_DIHEDRAL":
                dihedrals_top += "{:>6} {:>6} {:>6} {:>6}{:>2}{:>18.9e}{:>18.9e}{:>3d}\n".format(
                                dih.atmi.index + 1, dih.atmj.index + 1,
                                dih.atmk.index + 1, dih.atml.index + 1,
                                func, dih.phi0, dih.kd, dih.mult)
            elif dih.prefix_label == "HARMONIC_DIHEDRAL":
                dihedrals_top += "{:>6} {:>6} {:>6} {:>6}{:>2}{:>18.9e}{:>18.9e}\n".format(
                                dih.atmi.index + 1, dih.atmj.index + 1,
                                dih.atmk.index + 1, dih.atml.index + 1,
                                func, dih.phi0, dih.kd)
            else:
                print "Warning: unknown dihedral interaction for: {}".format(dih.describe())
        return dihedrals_top

    def _get_pairs_top(self):
        """ Get the [ pairs ] top for SBM Gromacs """
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
                                pot.eps, pot.r0, pot.width, pot.rNC)
                else:
                    print "Warning: interaction is not supported: {}".format(pot.describe())
                pairs_top += "{}{:>2}{}\n".format(atm_idxs, func, params)
        return pairs_top
    
    def _get_exclusions_top(self):
        """ Get [ exclusions ] top""" 
        exclusions_top = " [ exclusions ]\n"
        for pot in self.model.Hamiltonian.pairs:
            exclusions_top += "{:>6} {:>6}\n".format(pot.atmi.index + 1, pot.atmj.index + 1)
        return exclusions_top

    def generate_topology(self):
        """ Return a structure-based topology file. Currently only for one molecule. """

        self._generate_interaction_tables()

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

