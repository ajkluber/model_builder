import numpy as np

import mdtraj as md

import util

from hamiltonian import Hamiltonian

class StructureBasedHamiltonian(Hamiltonian):

    def __init__(self):
        Hamiltonian.__init__(self)

    def describe(self):
        pass

    def add_sbm_backbone(self, Model):

        if not hasattr(Model,"ref_traj"):
            raise AttributeError("Need to set reference structure model.set_reference()")

        # How should user defined parameters be set? 
        if self._default_parameters == {}:
            util.default_sbm_parameters_warning()
            self._use_sbm_default_parameters() 

        if self._default_potentials == {}:
            util.default_sbm_potentials_warning()
            self._use_sbm_default_potentials() 

        self._add_sbm_bonds(Model)
        self._add_sbm_angles(Model)
        self._add_sbm_dihedrals(Model)

    def add_sbm_contacts(self, Model):
        self._add_sbm_contacts(Model)

    def add_sbm_potentials(self, Model):

        self.add_sbm_backbone(Model)
        self.add_sbm_contacts(Model)

    def _add_sbm_bonds(self, Model):
        top = Model.mapping.top

        if hasattr(Model,"ref_traj"):
            kb = self._default_parameters["kb"]
            code = self._default_potentials["bond"]
            xyz = Model.ref_traj.xyz[0]
            for atm1, atm2 in top.bonds:
                r0 = np.linalg.norm(xyz[atm1.index,:] - xyz[atm2.index,:])
                self._add_bond(code, atm1, atm2, kb, r0)
        else:
            util.missing_reference_warning()

    def _add_sbm_angles(self, Model):
        structure = Model.mapping

        if hasattr(Model,"ref_traj"):
            ka = self._default_parameters["ka"]
            code = self._default_potentials["angle"]
            for atm1, atm2, atm3 in structure._angles:
                idxs = np.array([[atm1.index, atm2.index, atm3.index]])
                theta0 = md.compute_angles(Model.ref_traj, idxs)[0][0]
                self._add_angle(code, atm1, atm2, atm3, ka, theta0)
        else:
            util.missing_reference_warning()

    def _add_sbm_dihedrals(self, Model):
        structure = Model.mapping

        if hasattr(Model,"ref_traj"):
            kd = self._default_parameters["kd"]
            code = self._default_potentials["dihedral"]
            for atm1, atm2, atm3, atm4 in structure._dihedrals:
                idxs = np.array([[atm1.index, atm2.index, atm3.index, atm4.index]])
                phi0 = md.compute_dihedrals(Model.ref_traj, idxs)[0][0]

                #if temp_phi > 0:
                #    phi0 = 2.*np.pi - temp_phi
                #else:
                #    phi0 = -temp_phi

                self._add_dihedral(code, atm1, atm2, atm3, atm4, kd, phi0, 1)
                self._add_dihedral(code, atm1, atm2, atm3, atm4, 0.5*kd, phi0, 3)

            kd = self._default_parameters["ka"]
            code = self._default_potentials["improper_dihedral"]
            for atm1, atm2, atm3, atm4 in structure._improper_dihedrals:
                # How does mdtraj angle correspond to gromacs? 
                idxs = np.array([[atm1.index, atm2.index, atm3.index, atm4.index]])
                phi0 = md.compute_dihedrals(Model.ref_traj, idxs)[0][0]
                self._add_dihedral(code, atm1, atm2, atm3, atm4, kd, phi0)
        else:
            util.missing_reference_warning()

    def _add_sbm_contacts(self, Model):
        """Add structure-based model contacts"""
    
        structure = Model.mapping
        if hasattr(Model,"ref_traj"):
            eps = self._default_parameters["eps"]
            code = self._default_potentials["contact"]
            xyz = Model.ref_traj.xyz[0]
            for atm1, atm2 in structure._contact_pairs:
                r0 = np.linalg.norm(xyz[atm1.index,:] - xyz[atm2.index,:])
                self._add_pair(code, atm1, atm2, eps, r0)
        else:
            util.missing_reference_warning()

    def _use_sbm_default_parameters(self):
        self._default_parameters = {"kb":20000., # kJ/(mol nm^2)
                                    "ka":40.,  # kJ/(mol deg^2)
                                    "kd":1.,     # kJ/mol
                                    "eps":1}     # kJ/(mol nm)

    def _use_sbm_default_potentials(self):
        self._default_potentials = {"bond":"HARMONIC_BOND",
                                "angle":"HARMONIC_ANGLE",
                                "dihedral":"COSINE_DIHEDRAL",
                                "improper_dihedral":"HARMONIC_DIHEDRAL",
                                "contact":"LJ1210"}

    def calc_native_nonative_pair_energy(self, traj, n_native_pairs, sum=True):
        """Energy for pair interactions

        Parameters
        ----------
        traj : mdtraj.Trajectory
        
        sum : bool (opt.)
            If sum=True return the total energy.
        """
        r = md.compute_distances(traj, self._pair_idxs)
        if sum:
            Enat = np.zeros(traj.n_frames, float)
            Enon = np.zeros(traj.n_frames, float)
        else:
            Enat = np.zeros((traj.n_frames, n_native_pairs), float)
            Enon = np.zeros((traj.n_frames, self.n_pairs - n_native_pairs), float)

        for i in range(n_native_pairs):
            if sum:
                Enat += self._pairs[i].V(r[:,i])
            else:
                Enat[:,i] = self._pairs[i].V(r[:,i])

        for i in range(n_native_pairs, self.n_pairs):
            if sum:
                Enon += self._pairs[i].V(r[:,i])
            else:
                Enon[:,i - n_native_pairs] = self._pairs[i].V(r[:,i])
        return Enat, Enon

