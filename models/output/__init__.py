
SUPPORTED_VERSIONS = ["4.5.4","4.6.5","4.6.5_sbm"]



class GromacsFiles(object):


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

        self.supported_pair_potentials = ["LJ1210", 
                                        "GAUSSIAN",
                                        "LJ12GAUSSIAN"]

    def _generate_interaction_tables(self):
        """Generates tables of user-defined potentials"""

        # Determine which interactions need to be tabled 
        for pot in self.model.Hamiltonian.pairs:
            if pot.prefix_label not in supported_pair_interactions:  

        self.tablep = self._get_LJ1210_table()
        self.LJtable = self.tablep
        r = np.arange(0, 20.0, 0.002)
        self.tables = []
        self.tablenames = []
        for i in range(self.n_tables):
            pair_indx = tabled_pairs[i]
            table_name = "table_b%d.xvg" % (i+1)
            self.tablenames.append(table_name)

            table = np.zeros((len(r),3),float)
            table[1:,0] = r[1:]
            table[10:,1] = self.pair_V[pair_indx](r[10:]) 
            table[10:,2] = -1.*self.pair_dV[pair_indx](r[10:]) 
            self.tables.append(table)

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

    def _get_atoms_top(self):
        """ Generate the [ atoms ] top."""
        top = self.model.mapping.topology
        atoms_top = " [ atoms ]\n"
        atoms_top += " ;nr  type  resnr residue atom  cgnr charge  mass\n"
        for atom in top.atoms:
            atoms_top += " {:>5d}{:>4}{:>8d}{:>5}{:>4}{:>8}{:>8.3f}{:>8.3f}\n".format(
                                atom.index + 1, atom.name, atom.residue.index + 1, 
                                atom.residue.name, atom.name, atom.index + 1, 0.0, 1.0)
        return atoms_top

    def _get_bonds_top(self):
        """ Generate the [ bonds ] top."""
        bonds_top = " [ bonds ]\n"
        bonds_top += " ; ai aj func r0(nm) Kb\n"
        for bond in self.model.Hamiltonian.bonds:
            func = self._bond_funcs[bond.prefix_label]
            bonds_top += "{:>6} {:>6}{:>2}{:>18.9e}{:>18.9e}\n".format(
                                bond.atmi.index + 1, bond.atmj.index + 1, func, bond.r0, bond.kb) 
        return bonds_top

    def _get_angles_top(self):
        """ Generate the [ angles ] top."""
        angles_top = " [ angles ]\n"
        angles_top += " ; ai  aj  ak  func  th0(deg)   Ka\n"
        for angle in self.model.Hamiltonian.angles:
            func = self._angle_funcs[angle.prefix_label]
            angles_top += "{:>6} {:>6} {:>6}{:>2}{:>18.9e}{:>18.9e}\n".format(
                            angle.atmi.index + 1, angle.atmj.index + 1, 
                            angle.atmk.index + 1, func, 
                            angle.theta0, angle.ka)

    def _get_dihedrals_top(self):
        """ Generate the [ dihedrals ] top."""
        dihedrals_top = " [ dihedrals ]\n"
        dihedrals_top += " ; ai  aj  ak al  func  phi0(deg)   Kd mult\n"
        for dih in self.model.Hamiltonian.dihedrals:
            func = self._dihedral_funcs[dih.prefix_label]
            #self._dihedrals_top += "%6d %6d %6d %6d%2d%18.9e%18.9e%2d\n" %  \
            if dih.prefix_label == "COSINE_DIHEDRAL":
                dihedrals_top += "{:>6} {:>6} {:>6} {:>6}{:>2}{:>18.9e}{:>18.9e}\n".format(
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
        """ Get the [ pairs ] string for SBM Gromacs """
        pairs_string = " [ pairs ]\n"
        pairs_string += " ; %5s  %5s %4s    %8s   %8s\n" % ("i","j","type","c10","c12")
        for i in range(self.n_pairs):
            if i not in self.smog_pair_indxs:
                res_a = self.pairs[i][0]
                res_b = self.pairs[i][1]
                r0 = self.pairwise_other_parameters[i][0]
                eps = self.pair_eps[i]
                if self.pairwise_type[i] == 1:     # LJ12
                    c12 = eps*(r0**12)
                    c10 = 0
                    pairs_string += "%6d %6d%2d%18.9e%18.9e\n" % (res_a,res_b,1,c10,c12)
                elif self.pairwise_type[i] == 2:   # LJ1210
                    c12 = eps*5.0*(r0**12)
                    c10 = eps*6.0*(r0**10)
                    pairs_string += "%6d %6d%2d%18.9e%18.9e\n" % (res_a,res_b,1,c10,c12)
                else:
                    pass

    def _get_smog_pairs_top(self):
        pairs_string = " [ pairs ]\n"
        pairs_string += " ; Using smog Guassian interactions\n"
        for i in range(len(self.smog_pairs)):
            res_a = self.smog_pairs[i][0]
            res_b = self.smog_pairs[i][1]
            pairs_string += "%5d%5d%s\n" % (res_a,res_b,self.smog_strings[i]) 

    def generate_topology(self):
        """ Return a structure-based topology file. Currently only for one molecule. """

        top =  " ; Structure-based  topology file for Gromacs:\n"
        top += " [ defaults ]\n"
        top += " ;nbfunc comb-rule gen-pairs\n"
        top += "      1           1 no\n\n"
        #top += self.atomtypes_top #TODO
        top += " [ moleculetype ]\n"
        top += " ;name   nrexcl\n"
        top += " Macromolecule           3\n\n"

        top += "{}\n".format(self._get_atoms_top())
        top += "{}\n".format(self._get_bonds_top())
        #top += "{}\n".format(self._get_tabled_top()) #TODO
        top += "{}\n".format(self._get_angles_top())
        top += "{}\n".format(self._get_dihedrals_top())
        #top += "{}\n".format(self._get_pairs_top()) # TODO
        #top += "{}\n".format(self._get_exclusions_top()) # TODO

        top += " [ system ]\n"
        top += " ; name\n"
        top += " Macromolecule\n\n"
        top += " [ molecules ]\n"
        top += " ; name molec \n"
        top += " Macromolecule 1\n\n"

        self.topfile = top
