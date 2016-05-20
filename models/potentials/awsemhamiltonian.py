import numpy as np

import mdtraj as md

#from hamiltonian import Hamiltonian
from model_builder.models.mappings import AwsemBackboneMapping

import util
import awsem

class AwsemHamiltonian(object):

    def __init__(self):
        self.gamma_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
                               'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                               'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
                               'SER', 'THR', 'TRP', 'TYR', 'VAL']

        self.potential_types = [ "DIRECT", "WATER", "BURIAL" ]
        self.potential_gen = [ awsem.AWSEM_POTENTIALS[x] for x in self.potential_types ]
        self.potential_forms = { x:None for x in self.potential_types }

    @property
    def top(self):  
        return self.topology

    def _potential_is_parameterized(self, code):
        if self.potential_forms[code] is None:
            #print "Potential {} is not parameterized!".format(code)
            status = False
        else:
            status = True
        return True

    def _source_parameters(self, param_path):
        self._param_path = param_path 
        #self._backbone_coeff_path = backbone_coeff_path
        self._source_gammas()
        self._source_backbone_coeff()

        # Parameterize functional forms

    def _source_backbone_coeff(self):
        """Extract parameters from file
        
        AWSEM grabs parameters from a file called fix_backbone_coeff.data which
        control non-residue-specific properties of the force field, such as the
        boundaries of the contact well, the equilibrium distances between atoms
        in a hydrogen bond, etc.

        Almost all of these parameters are should stay at their default values.

        Distance units need to be converted from Ang. to nm.
        """
        with open("{}/fix_backbone_coeff.data".format(self._param_path), "r") as fin: 
            all_lines = fin.readlines()
            for i in range(len(all_lines)):
                line = all_lines[i]
                if line == "[Water]\n":
                    self._source_water_params(all_lines[i + 1:i + 8])
                elif line == "[Burial]\n":
                    self._source_burial_params(all_lines[i + 1:i + 6])

                # TODO:
                # - source helix params
                # - source rama params
                # - source dssp params
                # - source debye params

    def _source_water_params(self, lines):
        """Parameterize the form of the direct and water interaction"""
        # Convert distance units from Ang. to nm.  
        lambda_water = lambda_direct = float(lines[0])
        nu, nu_sigma = [ 10.*float(x) for x in lines[1].split() ]
        direct_r_min, direct_r_max, dum = [ float(x)/10. for x in lines[5].split() ]
        water_r_min, water_r_max, dum = [ float(x)/10. for x in lines[6].split() ]

        self.potential_forms["DIRECT"] = awsem.AWSEM_POTENTIALS["DIRECT"](
                lambda_direct=lambda_direct, nu=nu, 
                r_min=direct_r_min, r_max=direct_r_max)

        self.potential_forms["WATER"] = awsem.AWSEM_POTENTIALS["WATER"](
                lambda_water=lambda_water, nu=nu, nu_sigma=nu_sigma, 
                r_min=water_r_min, r_max=water_r_max)

    def _source_burial_params(self, lines):
        """Parameterize the form of the burial interaction"""
        # parameterize the burial interaction.                    
        lambda_burial = float(lines[0])
        nu = float(lines[1])
        rho1_lims = [ float(x) for x in lines[2].split() ]
        rho2_lims = [ float(x) for x in lines[3].split() ]
        rho3_lims = [ float(x) for x in lines[4].split() ]

        self.potential_forms["BURIAL"] = awsem.AWSEM_POTENTIALS["BURIAL"](
                lambda_burial=lambda_burial, nu=nu, rho1_lims=rho1_lims, 
                rho2_lims=rho2_lims, rho3_lims=rho3_lims)

    def _source_gammas(self):
        """Extract contact and burial gamma parameters from file"""

        self.gamma_burial = np.loadtxt("{}/burial_gamma.dat".format(self._param_path))

        with open("{}/gamma.dat".format(self._param_path), "r") as fin:
            gamma_direct = np.zeros((20,20))
            for i in range(20):
                for j in range(i, 20):
                    line = fin.readline()
                    gamma_direct[i, j] = gamma_direct[j, i] = float(line.split()[0])

            line = fin.readline()
            gamma_water = np.zeros((20,20))
            gamma_protein = np.zeros((20,20))
            for i in range(20):
                for j in range(i, 20):
                    line = fin.readline()
                    gamma_water[i, j] = gamma_water[j,i] = float(line.split()[0])
                    gamma_protein[i, j] = gamma_protein[j, i] = float(line.split()[1])

        self.gamma_direct = gamma_direct 
        self.gamma_water = gamma_water 
        self.gamma_protein = gamma_protein 

    def set_topology(self, topology, parameterize=False):
        """Set the topology used to construct Hamiltonian parameters
        
        Parameters
        ----------
        topology : mdtraj.Topology
            Coarse-grained topology corresponding to a AwsemMapping.
        """
        self.backbone_mapping = AwsemBackboneMapping(topology)
        self.topology = self.backbone_mapping.topology

        if parameterize:
            # Set the contact, burial potentials using this topology.
            self.parameterize()

    def parameterize(self):
        """Parameterize the potential terms based off the current topology
        
        Extract the residue identites and atomic indices from the topology in
        order to assign the transferable parameters.
        """

        self._parameterize_burial()
        self._parameterize_contacts()
        #self._parameterize_alpha_helical()

    def _parameterize_burial(self):
        res_burial = []
        for i in range(self.top.n_residues):
            # assign burial parameters
            resi = self.top.residue(i)
            res_idx1 = self.gamma_residues.index(resi.name)
            res_burial.append(self.gamma_burial[res_idx1])
        self.res_gamma_burial = np.array(res_burial)

    def _parameterize_contacts(self):
        """Assign gammas for each contact interaction"""

        # Pre-allocate?
        n_pairs = self.top.n_residues*(self.top.n_residues - 1)
        contact_pairs = []
        contact_gamma_direct = []
        contact_gamma_water = []
        contact_gamma_protein = []
        for i in range(self.top.n_residues):
            resi = self.top.residue(i)
            res_idx1 = self.gamma_residues.index(resi.name)
            if resi.name == "GLY":
                idx1 = self.top.select("resid {} and name CA".format(i))[0]
            else:
                idx1 = self.top.select("resid {} and name CB".format(i))[0]
            for j in range(self.top.n_residues):
                if i != j:
                    resj = self.top.residue(j)

                    if resj.name == "GLY":
                        idx2 = self.top.select("resid {} and name CA".format(j))[0]
                    else:
                        idx2 = self.top.select("resid {} and name CB".format(j))[0]

                    # track the atom indices involved in the contact.
                    contact_pairs.append([idx1, idx2])

                    res_idx2 = self.gamma_residues.index(resj.name)

                    # Assign contact parameters for residues (i, j)
                    contact_gamma_direct.append(self.gamma_direct[res_idx1, res_idx2])
                    contact_gamma_water.append(self.gamma_water[res_idx1, res_idx2])
                    contact_gamma_protein.append(self.gamma_protein[res_idx1, res_idx2])

        self.contact_pairs = np.array(contact_pairs)
        self.contact_gamma_direct = np.array(contact_gamma_direct)
        self.contact_gamma_water = np.array(contact_gamma_water)
        self.contact_gamma_protein = np.array(contact_gamma_protein)
        self.n_pairs = len(self.contact_pairs) 

    def _parameterize_alpha_helical(self):
        # indices for alpha helical hydrogen bonding terms and their strengths
        pass


    def calculate_direct_energy(self, traj, sum=True):
        """Calculate the two-body direct contact potential
        
        Parameters
        ----------
        traj : mdtraj.Trajectory
            Trajectory to calculate energy over.
        sum : opt, bool
            If true (default) return the sum of the burial potentials. If
            false, return the burial energy of each individual residue.
        """

        bb_traj = self.backbone_mapping.map_traj(traj)
        direct = self.potential_forms["DIRECT"] 
        r = md.compute_distances(bb_traj, self.contact_pairs)

        if sum:
            Vdirect = np.zeros(bb_traj.n_frames, float)
        else:
            Vdirect = np.zeros((bb_traj.n_frames, self.n_pairs), float)

        for i in range(self.n_pairs):
            if sum:
                Vdirect += direct.V(r[:,i], self.contact_gamma_direct[i])
            else:
                Vdirect[:,i] = direct.V(r[:,i], self.contact_gamma_direct[i])
        return Vdirect 

    def calculate_water_energy(self, traj, local_density=None, sum=True):
        """Calculate the one-body burial potential
        
        Parameters
        ----------
        traj : mdtraj.Trajectory
            Trajectory to calculate energy over.
        sum : opt, bool
            If true (default) return the sum of the burial potentials. If
            false, return the burial energy of each individual residue.
        """
        pass

    def calculate_burial_energy(self, traj, local_density=None, sum=True):
        """Calculate the one-body burial potential
        
        Parameters
        ----------
        traj : mdtraj.Trajectory
            Trajectory to calculate energy over.
        local_density : np.ndarray (traj.n_frames, traj.n_residues)
            Local protein density around each residue for all frames in traj. 
        sum : opt, bool
            If true (default) return the sum of the burial potentials. If
            false, return the burial energy of each individual residue.
        """

        if local_density is None:
            bb_traj = self.backbone_mapping.map_traj(traj)
            res_local_density = self._calculate_local_density(bb_traj)
        else:
            res_local_density = local_density

        if sum:
            Vburial = np.zeros(traj.n_frames, float)
        else:
            Vburial = np.zeros((traj.n_frames, self.top.n_residues), float)

        burial = self.potential_forms["BURIAL"]
        for i in range(self.top.n_residues):
            if sum:
                Vburial += burial.V(res_local_density[:,i], self.res_gamma_burial[i])
            else:
                Vburial[:,i] = burial.V(res_local_density[:,i], self.res_gamma_burial[i])

        return Vburial

    def _calculate_local_density(self, traj):
        """Calculate local protein density around each residue"""

        direct = self.potential_forms["DIRECT"]
        contact_dists = md.compute_distances(traj, self.contact_pairs)
        contacts = direct.theta_I(contact_dists)
        contact_matrix = contacts.reshape((traj.n_frames, traj.n_residues, traj.n_residues - 1))
        local_density = np.sum(contact_matrix, axis=2)

        return local_density

    def calculate_energy(self):
        pass

if __name__ == "__main__":
    pass
