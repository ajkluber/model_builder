import numpy as np

import potentials

from model import Model

class StructureBasedModel(Model):

    def __init__(self, topology, bead_repr=None):
        """Structure-based Model (SBM)

        Parameters
        ----------
        topology : mdtraj.Topology object
            An mdtraj Topology object that describes the molecular topology.

        bead_repr : str [CA, CACB]
            A code specifying the desired coarse-grain mapping. The all-atom
        to coarse-grain mapping.

        Methods
        -------
        assign_* :
            Methods assign which atoms have bonded constraints
        (angle potentials, dihedral, etc.)

        add_* :
            Methods add potentials to the Hamiltonian.

        """



        Model.__init__(self, topology, bead_repr=bead_repr)
        self.Hamiltonian = potentials.StructureBasedHamiltonian()
        self.mapping.add_atoms()

    def set_reference(self, traj):
        """Set the reference structure

        Parameters
        ----------
        traj : mdtraj.Trajectory object
            Trajectory to be used as a reference (only uses first frame). The
            geometry of the reference structure is used to construct structure-
            based potentials.
        """
        self.ref_traj_aa = traj[0]
        self.ref_traj = self.mapping.map_traj(traj[0])

    def save_starting_conf(self, saveas="conf.gro"):
        if hasattr(self, "starting_traj"):
            self.starting_traj.save(saveas)
        else:
            self.ref_traj.save(saveas)

    def assign_disulfides(self, disulfides, simple=False):
        self.mapping.add_disulfides(disulfides, simple=simple)

    def assign_backbone(self):
        self.mapping._assign_sbm_angles()
        self.mapping._assign_sbm_dihedrals()

    def assign_contacts(self):
        self.mapping._assign_sbm_contacts(self.ref_traj_aa)

    def add_sbm_potentials(self):
        self.Hamiltonian.add_sbm_potentials(self)

    def add_sbm_backbone(self):
        self.Hamiltonian.add_sbm_backbone(self)

    def add_sbm_contacts(self):
        self.Hamiltonian._add_sbm_contacts(self)

    def use_flavored_nonnative_interactions(self, type_dict):
        self.mapping._assign_flavored_nonnative_values(type_dict)
