
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
