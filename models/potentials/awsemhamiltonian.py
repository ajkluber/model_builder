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
                if line == "[Burial]\n":
                    self._source_burial_params(all_lines[i + 1:i + 6])
                elif line == "[Water]\n":
                    self._source_water_params(all_lines[i + 1:i + 8])
                elif line == "[Helix]\n":
                    self._source_helix_params(all_lines[i + 1:i + 12])

                # TODO:
                # - source rama params
                # - source dssp params
                # - source debye params

    def _source_burial_params(self, lines):
        """Parameterize the form of the burial interaction"""

        lambda_burial = float(lines[0])
        nu = float(lines[1])
        rho1_lims = [ float(x) for x in lines[2].split() ]
        rho2_lims = [ float(x) for x in lines[3].split() ]
        rho3_lims = [ float(x) for x in lines[4].split() ]

        self.potential_forms["BURIAL"] = awsem.AWSEM_POTENTIALS["BURIAL"](
                lambda_burial=lambda_burial, nu=nu, rho1_lims=rho1_lims, 
                rho2_lims=rho2_lims, rho3_lims=rho3_lims)

    def _source_water_params(self, lines):
        """Parameterize the form of the direct and water interaction"""

        lambda_water = lambda_direct = float(lines[0])
        nu, nu_sigma = [ float(x) for x in lines[1].split() ]
        nu *= 10.
        rho_0 = float(lines[2])
        self.contact_exclude_neighbors = int(lines[3])
        direct_r_min, direct_r_max, dum = [ float(x)/10. for x in lines[5].split() ]
        water_r_min, water_r_max, dum = [ float(x)/10. for x in lines[6].split() ]

        self.potential_forms["DIRECT"] = awsem.AWSEM_POTENTIALS["DIRECT"](
                lambda_direct=lambda_direct, nu=nu, 
                r_min=direct_r_min, r_max=direct_r_max)

        self.potential_forms["WATER"] = awsem.AWSEM_POTENTIALS["WATER"](
                lambda_water=lambda_water, nu=nu, nu_sigma=nu_sigma, 
                r_min=water_r_min, r_max=water_r_max, rho_0=rho_0)

    def _source_helix_params(self, lines):
        """Parameterize the form of the helix interaction"""

        lambda_helix = float(lines[0])
        gamma_protein, gamma_water = [ float(x) for x in lines[1].split() ]
        nu, nu_sigma = [ float(x) for x in lines[2].split() ]
        nu *= 10.
        rho_0 = float(lines[3]) 
        helix_i_diff = int(lines[4]) 
        helix_cutoff = float(lines[5])/10.
        self.res_helix_fai = [ float(x) for x in lines[7].split() ]
        pro_vals = [ x for x in lines[8].split() ]
        self._pro_acceptor_flag = int(pro_vals[0])
        self._pro_acceptor_fai = float(pro_vals[0]) 
        sigma_OH, sigma_ON = [ float(x)/10. for x in lines[9].split() ]
        r_OH, r_ON = [ float(x)/10. for x in lines[10].split() ]

        self.potential_forms["HELIX"] = awsem.AWSEM_POTENTIALS["HELIX"](
                lambda_helix=lambda_helix, gamma_protein=gamma_protein,
                gamma_water=gamma_water, nu=nu, nu_sigma=nu_sigma, rho_0=rho_0,
                r_ON=r_ON, r_OH=r_OH, sigma_ON=sigma_ON, sigma_OH=sigma_OH)

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
                    gamma_water[i, j] = gamma_water[j, i] = float(line.split()[1])
                    gamma_protein[i, j] = gamma_protein[j, i] = float(line.split()[0])

        self.gamma_direct = gamma_direct 
        self.gamma_water = gamma_water 
        self.gamma_protein = gamma_protein 

    def set_topology(self, topology, parameterize=False):
        """Set the topology used to construct Hamiltonian parameters

        
        Parameters
        ----------
        topology : mdtraj.Topology
            Coarse-grained topology corresponding to an AwsemMapping
            representation i.e. three atoms per residue. This is mapped to a
            representation that includes the backbone atoms
            (AwsemBackboneMapping), so that hydrogen bonding can be calculated.
            The indices used for calculating all energy terms are in terms of
            the full backbone topology.
        """
        self.backbone_mapping = AwsemBackboneMapping(topology)
        self.topology = self.backbone_mapping.topology

        if parameterize:
            # Set the contact, burial, helix potentials using this topology.
            self.parameterize()

    def parameterize(self):
        """Parameterize the potential terms based off the current topology
        
        Extract the residue identites and atomic indices from the topology in
        order to assign the transferable parameters.
        """

        self._parameterize_burial()
        self._parameterize_contacts()
        self._parameterize_alpha_helical()

    def _parameterize_burial(self):
        res_burial = []
        for i in range(self.top.n_residues):
            # assign burial parameters
            resi = self.top.residue(i)
            res_idx1 = self.gamma_residues.index(resi.name)
            res_burial.append(self.gamma_burial[res_idx1])
        self.res_gamma_burial = np.array(res_burial)

        # pairs included in the local density calculation for each residue. We
        # exclude -/+1 neighbors on the same chain.
        self.res_burial_pairs = []
        for i in range(self.top.n_residues):
            temp_res_burial_pairs = []
            resi = self.top.residue(i)
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

                    if (resi.chain.index == resj.chain.index):
                        if abs(resj.index - resi.index) > 1:
                            temp_res_burial_pairs.append([idx1, idx2])
                    else:
                        temp_res_burial_pairs.append([idx1, idx2])

            self.res_burial_pairs.append(np.array(temp_res_burial_pairs))

    def _parameterize_contacts(self):
        """Assign gammas for each contact interaction"""

        contact_pairs = []
        contact_res_idxs = []
        contact_gamma_idxs = []
        for i in range(self.top.n_residues):
            resi = self.top.residue(i)
            gamma_idx1 = self.gamma_residues.index(resi.name)

            if resi.name == "GLY":
                idx1 = self.top.select("resid {} and name CA".format(i))[0]
            else:
                idx1 = self.top.select("resid {} and name CB".format(i))[0]
            for j in range(i + 1, self.top.n_residues):
                resj = self.top.residue(j)
                gamma_idx2 = self.gamma_residues.index(resj.name)

                if resj.name == "GLY":
                    idx2 = self.top.select("resid {} and name CA".format(j))[0]
                else:
                    idx2 = self.top.select("resid {} and name CB".format(j))[0]

                if (resi.chain.index == resj.chain.index):
                    # exclude neighbors on same chain.
                    if (resj.index - resi.index) >= self.contact_exclude_neighbors:
                        contact_pairs.append([idx1, idx2])
                        contact_res_idxs.append([resi.index, resj.index])
                        contact_gamma_idxs.append([gamma_idx1, gamma_idx2])
                else:
                    contact_pairs.append([idx1, idx2])
                    contact_res_idxs.append([resi.index, resj.index])
                    contact_gamma_idxs.append([gamma_idx1, gamma_idx2])

        self.n_pairs = len(contact_pairs)
        self._contact_pairs = np.array(contact_pairs)
        self._contact_res_idxs = np.array(contact_res_idxs)
        self._contact_gamma_idxs = np.array(contact_gamma_idxs)

    def _parameterize_alpha_helical(self):
        # indices for alpha helical hydrogen bonding terms and their strengths

        helix_res_idxs = []
        helix_ON_pairs = []
        helix_OH_pairs = []
        helix_fai = []
        helix_fai_4 = []
        # indices participating in alpha-helical hydrogen bonds
        for i in range(self.top.n_residues - 4):
            resi = self.top.residue(i) 
            resi_4 = self.top.residue(i + 4) 
            if resi.chain.index == resi_4.chain.index:
                # Get donor acceptor indices
                O_idx = self.top.select("resid {} and name O".format(resi.index))[0]
                H_idx = self.top.select("resid {} and name H".format(resi_4.index))[0] 
                N_idx = self.top.select("resid {} and name N".format(resi_4.index))[0]
                helix_res_idxs.append([resi.index, resi_4.index])
                helix_ON_pairs.append([O_idx, N_idx])
                helix_OH_pairs.append([O_idx, H_idx])
                gamma_idx1 = self.gamma_residues.index(resi.name)
                gamma_idx2 = self.gamma_residues.index(resi_4.name)

                # properly treat prolines. They can't be H-bond donors.
                if (resi.name == "PRO") and (self._pro_acceptor_flag == 1):
                    helix_fai.append(self._pro_acceptor_fai)
                else:
                    helix_fai.append(self.res_helix_fai[gamma_idx1])
                helix_fai_4.append(self.res_helix_fai[gamma_idx2])

        self._helix_res_idxs = np.array(helix_res_idxs)
        self._helix_ON_pairs = np.array(helix_ON_pairs)
        self._helix_OH_pairs = np.array(helix_OH_pairs)
        self._helix_fai = np.array(helix_fai)
        self._helix_fai_4 = np.array(helix_fai_4)
        self.n_alpha_helix = len(helix_ON_pairs)

    def _calculate_local_density(self, traj):
        """Calculate local protein density around each residue"""

        local_density = np.zeros((traj.n_frames, traj.n_residues))
        direct = self.potential_forms["DIRECT"]

        for i in range(traj.n_residues):
            # compute local density of residue
            pairs = self.res_burial_pairs[i]
            local_density[:,i] = np.sum(direct.theta_I(md.compute_distances(traj, pairs)), axis=1)
        return local_density

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
        burial = self.potential_forms["BURIAL"]

        if local_density is None:
            bb_traj = self.backbone_mapping.map_traj(traj)
            res_local_density = self._calculate_local_density(bb_traj)
        else:
            res_local_density = local_density

        if sum:
            Vburial = np.zeros(traj.n_frames, float)
        else:
            Vburial = np.zeros((traj.n_frames, self.top.n_residues), float)

        for i in range(self.top.n_residues):
            if sum:
                Vburial += burial.V(res_local_density[:,i], self.res_gamma_burial[i])
            else:
                Vburial[:,i] = burial.V(res_local_density[:,i], self.res_gamma_burial[i])

        return Vburial

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
        direct = self.potential_forms["DIRECT"] 

        bb_traj = self.backbone_mapping.map_traj(traj)
        r = md.compute_distances(bb_traj, self._contact_pairs)

        if sum:
            Vdirect = np.zeros(bb_traj.n_frames, float)
        else:
            Vdirect = np.zeros((bb_traj.n_frames, self.n_pairs), float)
        
        for i in range(self.n_pairs):
            gamma_direct = self.gamma_direct[self._contact_gamma_idxs[i,0], self._contact_gamma_idxs[i,1]]
            if sum:
                Vdirect += direct.V(r[:,i], gamma_direct)
            else:
                Vdirect[:,i] = direct.V(r[:,i], gamma_direct)
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
        water = self.potential_forms["WATER"] 

        bb_traj = self.backbone_mapping.map_traj(traj)
        r = md.compute_distances(bb_traj, self._contact_pairs)

        if local_density is None:
            res_local_density = self._calculate_local_density(bb_traj)
        else:
            res_local_density = local_density

        if sum:
            Vwater = np.zeros(bb_traj.n_frames, float)
        else:
            Vwater = np.zeros((bb_traj.n_frames, self.n_pairs), float)
        
        for i in range(self.n_pairs):
            rhoi = res_local_density[:,self._contact_res_idxs[i,0]]
            rhoj = res_local_density[:,self._contact_res_idxs[i,1]]
            gamma_water = self.gamma_water[self._contact_gamma_idxs[i,0], self._contact_gamma_idxs[i,1]]
            gamma_protein = self.gamma_protein[self._contact_gamma_idxs[i,0], self._contact_gamma_idxs[i,1]]
            if sum:
                Vwater += water.V(r[:,i], rhoi, rhoj, gamma_water, gamma_protein)
            else:
                Vwater[:,i] = water.V(r[:,i], rhoi, rhoj, gamma_water, gamma_protein)
        return Vwater 

    def calculate_helix_energy(self, traj, local_density=None, sum=True):
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
        helix = self.potential_forms["HELIX"]

        bb_traj = self.backbone_mapping.map_traj(traj)

        r_ON = md.compute_distances(bb_traj, self._helix_ON_pairs)
        r_OH = md.compute_distances(bb_traj, self._helix_OH_pairs)

        if local_density is None:
            res_local_density = self._calculate_local_density(bb_traj)
        else:
            res_local_density = local_density

        if sum:
            Vhelix = np.zeros(bb_traj.n_frames, float)
        else:
            Vhelix = np.zeros((bb_traj.n_frames, self.n_alpha_helix), float)

        for i in range(self.n_alpha_helix):
            rhoi = res_local_density[:,self._helix_res_idxs[i,0]]
            rhoi_4 = res_local_density[:,self._helix_res_idxs[i,1]]
            fai = self._helix_fai[i]
            fai_4 = self._helix_fai_4[i]
            if sum:
                Vhelix += helix.V(r_ON[:,i], r_OH[:,i], rhoi, rhoi_4, fai, fai_4)
            else:
                Vhelix[:, i] = helix.V(r_ON[:,i], r_OH[:,i], rhoi, rhoi_4, fai, fai_4)
        return Vhelix

    def calculate_energy(self):
        pass

if __name__ == "__main__":
    pass
