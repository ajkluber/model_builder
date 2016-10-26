import numpy as np
import os

import mdtraj as md

#from hamiltonian import Hamiltonian
from model_builder.models.mappings import AwsemBackboneMapping

import util
import awsem

class AwsemHamiltonian(object):

    def __init__(self, three_bead_topology):
        """AWSEM model

        Parameters
        ----------
        three_bead_topology : mdtraj.Topology
            Topology of the coarse-grain AWSEM representation. AWSEM is
            simulated with three explicit atoms per residue (CA, CB, O) but the
            potentials depend on the three remaining backbone atoms (C, N, H)
            that can be geometrically interpolated. All of the indices used
            within this Hamiltonian are from the full-backbone topology.

        """

        self.gamma_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                               'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                               'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                               'SER', 'THR', 'TRP', 'TYR', 'VAL']

        self.potential_types = [ "DIRECT", "WATER", "BURIAL" ]
        self.potential_gen = [ awsem.AWSEM_POTENTIALS[x] for x in self.potential_types ]
        self.potential_forms = { x:None for x in self.potential_types }

        self.three_bead_topology = three_bead_topology
        self.backbone_mapping = AwsemBackboneMapping(three_bead_topology)
        self.topology = self.backbone_mapping.topology

        self._get_terminal_residues()

    @property
    def top(self):
        return self.topology

    def _set_charged_residues(self, charged_residues):
        self._charged_residues = charged_residues

    def _potential_is_parameterized(self, code):
        if self.potential_forms[code] is None:
            #print "Potential {} is not parameterized!".format(code)
            status = False
        else:
            status = True
        return status

    def _get_terminal_residues(self):
        """Get which residues are on the end of a chain"""
        self._n_terminal_residues = []
        self._c_terminal_residues = []
        for chain in self.top.chains:
            self._n_terminal_residues.append(chain.residue(0).index)
            self._c_terminal_residues.append(chain.residue(-1).index)

    def _source_parameters(self, param_path):
        self._param_path = param_path
        self._source_gammas()
        self._source_backbone_coeff()

        # Parameterize functional forms
        self._parameterize()

    def _source_backbone_coeff(self):
        """Extract parameters from file

        AWSEM grabs parameters from a file called fix_backbone_coeff.data which
        control non-residue-specific properties of the force field, such as the
        boundaries of the contact well, the equilibrium distances between atoms
        in a hydrogen bond, etc.

        Almost all of these parameters are should stay at their default values.

        Distance units need to be converted from Ang. to nm.
        """

        #set default parameters first, then override them with the
        self._set_default_backbone_coeff()

        with open("{}/fix_backbone_coeff.data".format(self._param_path), "r") as fin:
            all_lines = fin.readlines()
            for i in range(len(all_lines)):
                line = all_lines[i]
                if line == "[Burial]\n":
                    self._source_burial_params(all_lines[i + 1:i + 6])
                elif line == "[Water]\n":
                    self._source_water_params(all_lines[i + 1:i + 8])
                elif line == "[DebyeHuckel]\n":
                    self._source_debye_params(all_lines[i + 1:i + 5])
                elif line == "[Helix]\n":
                    self._source_helix_params(all_lines[i + 1:i + 12])
                elif line == "[Rama]\n":
                    self._source_rama_params(all_lines[i + 1:i + 8])
                elif line == "[Rama_P]\n":
                    self._source_rama_proline_params(all_lines[i + 1:i + 7])

                # TODO:
                # - source rama params
                # - source dssp params
                # - source debye params

    def _set_default_backbone_coeff(self):
        """ Set the default backbone parameters """

        #Burial Params
        self.potential_forms["BURIAL"] = awsem.AWSEM_POTENTIALS["BURIAL"](
                lambda_burial=1.0, nu=4.0, rho1_lims=[0.0, 3.0],
                rho2_lims=[3.0, 6.0], rho3_lims=[6.0, 9.0])

        #Water params
        self.potential_forms["DIRECT"] = awsem.AWSEM_POTENTIALS["DIRECT"](
                lambda_direct=10., nu=50., r_min=0.45, r_max=0.65)

        self.potential_forms["WATER"] = awsem.AWSEM_POTENTIALS["WATER"](
                lambda_water=10., nu=50., nu_sigma=7.0, r_min=0.65,
                r_max=0.95, rho_0=2.6)

        self._contact_exclude_neighbors = 13

        #DebyeHuckel params from amylometer branch on github
        self._debye_kplusplus = 0.1
        self._debye_kminusminus = 0.1
        self._debye_kplusminus = 0.1
        self._debye_exclude_neighbors = 10

        self.potential_forms["DEBYE"] = awsem.AWSEM_POTENTIALS["DEBYE"](
                k_screening=1.0, debye_length=1.0)

        #Helix params
        self.potential_forms["HELIX"] = awsem.AWSEM_POTENTIALS["HELIX"](
                lambda_helix=1.5, gamma_protein=2.0,
                gamma_water=-1.0, nu=70., nu_sigma=7.0, rho_0=3.0,
                r_ON=.298, r_OH=.206, sigma_ON=0.068, sigma_OH=0.076)
        self.res_helix_fai = [0.77, 0.68, 0.07, 0.15, 0.23, 0.33, 0.27,
                0.0, 0.06, 0.23, 0.62, 0.65, 0.50, 0.41, -3.0, 0.35,
                0.11, 0.45, 0.17, 0.14]
        self._pro_acceptor_flag = 0
        self._pro_acceptor_fai = 0.0

        #Rama params
        rama_fields = np.array([ [1.3149, 15.398, 0.15, 1.74, 0.65, -2.138],
                [1.32016, 49.0521, 0.25, 1.265, 0.45, 0.318],
                [1.0264, 49.0954, 0.65, -1.041, 0.25, -0.78] ])

        W = rama_fields[:, 0]
        sigma = rama_fields[:, 1]
        omega_phi = rama_fields[:, 2]
        phi0 = -rama_fields[:, 3]
        omega_psi = rama_fields[:, 4]
        psi0 = -rama_fields[:, 5]

        self.potential_forms["RAMA"] = awsem.AWSEM_POTENTIALS["RAMA"](
                lambda_rama=2.0, W=W, sigma=sigma,
                omega_phi=omega_phi, phi0=phi0,
                omega_psi=omega_psi, psi0=psi0)

        alpha = [2.0, 419.0, 1.0, 0.995, 1.0, 0.820]
        beta = [2.0, 15.398, 1.0, 2.25, 1.0, -2.16]

        self.potential_forms["RAMA_ALPHA"] = awsem.AWSEM_POTENTIALS["RAMA"](
                lambda_rama=2.0, W=alpha[0], sigma=alpha[1],
                omega_phi=alpha[2], phi0=-alpha[3],
                omega_psi=alpha[4], psi0=-alpha[5])

        self.potential_forms["RAMA_BETA"] = awsem.AWSEM_POTENTIALS["RAMA"](
                lambda_rama=2.0, W=beta[0], sigma=beta[1],
                omega_phi=beta[2], phi0=-beta[3],
                omega_psi=beta[4], psi0=-beta[5])

        #Rama_P params
        rama_fields = np.array([ [0.0, 0.0, 1.0, 0.0, 1.0, 0.0],
                [2.17, 105.52, 1.0, 1.153, 0.15, -2.4],
                [2.15, 109.09, 1.0, 0.95, 0.15, 0.218] ])
        self.potential_forms["RAMA_PROLINE"] = awsem.AWSEM_POTENTIALS["RAMA"](
                W=W, sigma=sigma, omega_phi=omega_phi, phi0=phi0,
                omega_psi=omega_psi, psi0=psi0)


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

    def _source_water_params(self, lines):
        """Parameterize the form of the direct and water interaction"""

        lambda_water = lambda_direct = float(lines[0])
        nu, nu_sigma = [ float(x) for x in lines[1].split() ]
        nu *= 10.
        rho_0 = float(lines[2])
        self._contact_exclude_neighbors = int(lines[3])
        direct_r_min, direct_r_max, dum = [ float(x)/10. for x in lines[5].split() ]
        water_r_min, water_r_max, dum = [ float(x)/10. for x in lines[6].split() ]

        self.potential_forms["DIRECT"] = awsem.AWSEM_POTENTIALS["DIRECT"](
                lambda_direct=lambda_direct, nu=nu,
                r_min=direct_r_min, r_max=direct_r_max)

        self.potential_forms["WATER"] = awsem.AWSEM_POTENTIALS["WATER"](
                lambda_water=lambda_water, nu=nu, nu_sigma=nu_sigma,
                r_min=water_r_min, r_max=water_r_max, rho_0=rho_0)

    def _source_debye_params(self, lines):

        self._debye_kplusplus, self._debye_kminusminus, self._debye_kplusminus =\
                [ float(x)/10. for x in lines[0].split() ]
        k_screening = float(lines[1])
        debye_length = float(lines[2])/10.
        self._debye_exclude_neighbors = int(lines[3])

        self.potential_forms["DEBYE"] = awsem.AWSEM_POTENTIALS["DEBYE"](
                k_screening=k_screening, debye_length=debye_length)

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

    def _source_rama_params(self, lines):
        """Parameterize the form of the rama interaction"""

        lambda_rama = float(lines[0])
        n_rama_fields = int(lines[1])
        rama_fields = np.array([ [ float(x) for x in line_temp.split() ] for line_temp in lines[2:5] ])
        W = rama_fields[:, 0]
        sigma = rama_fields[:, 1]
        omega_phi = rama_fields[:, 2]
        phi0 = -rama_fields[:, 3]
        omega_psi = rama_fields[:, 4]
        psi0 = -rama_fields[:, 5]

        self.potential_forms["RAMA"] = awsem.AWSEM_POTENTIALS["RAMA"](
                lambda_rama=lambda_rama, W=W, sigma=sigma,
                omega_phi=omega_phi, phi0=phi0,
                omega_psi=omega_psi, psi0=psi0)

        # parameterize secondary structure bias terms
        alpha = [ float(x) for x in lines[5].split() ]
        beta = [ float(x) for x in lines[6].split() ]

        self.potential_forms["RAMA_ALPHA"] = awsem.AWSEM_POTENTIALS["RAMA"](
                lambda_rama=lambda_rama, W=alpha[0], sigma=alpha[1],
                omega_phi=alpha[2], phi0=-alpha[3],
                omega_psi=alpha[4], psi0=-alpha[5])

        self.potential_forms["RAMA_BETA"] = awsem.AWSEM_POTENTIALS["RAMA"](
                lambda_rama=lambda_rama, W=beta[0], sigma=beta[1],
                omega_phi=beta[2], phi0=-beta[3],
                omega_psi=beta[4], psi0=-beta[5])

    def _source_rama_proline_params(self, lines):
        """Parameterize the form of the rama interaction"""

        n_rama_fields = int(lines[0])
        rama_fields = np.array([ [ float(x) for x in line_temp.split() ] for line_temp in lines[1:4] ])
        W = rama_fields[:, 0]
        sigma = rama_fields[:, 1]
        omega_phi = rama_fields[:, 2]
        phi0 = -rama_fields[:, 3]
        omega_psi = rama_fields[:, 4]
        psi0 = -rama_fields[:, 5]

        self.potential_forms["RAMA_PROLINE"] = awsem.AWSEM_POTENTIALS["RAMA"](
                W=W, sigma=sigma, omega_phi=omega_phi, phi0=phi0,
                omega_psi=omega_psi, psi0=psi0)

#    def _default_parameters(self):
#        """Create potential forms with default parameters"""
#
#        self.potential_forms["BURIAL"] = awsem.AWSEM_POTENTIALS["BURIAL"]()
#        self.potential_forms["DIRECT"] = awsem.AWSEM_POTENTIALS["DIRECT"]()
#        self.potential_forms["WATER"] = awsem.AWSEM_POTENTIALS["WATER"]()
#        self.potential_forms["HELIX"] = awsem.AWSEM_POTENTIALS["HELIX"]()
#        self.potential_forms["RAMA"] = awsem.AWSEM_POTENTIALS["RAMA"]()
#
#        self.potential_forms["RAMA_ALPHA"] = awsem.AWSEM_POTENTIALS["RAMA"](
#                lambda_rama=lambda_rama, W=alpha[0], sigma=alpha[1],
#                omega_phi=alpha[2], phi0=-alpha[3],
#                omega_psi=alpha[4], psi0=-alpha[5])
#
#        self.potential_forms["RAMA_BETA"] = awsem.AWSEM_POTENTIALS["RAMA"](
#                lambda_rama=lambda_rama, W=beta[0], sigma=beta[1],
#                omega_phi=beta[2], phi0=-beta[3],
#                omega_psi=beta[4], psi0=-beta[5])
#
#        self.potential_forms["RAMA_PROLINE"] = awsem.AWSEM_POTENTIALS["RAMA"](
#                W=W, sigma=sigma, omega_phi=omega_phi, phi0=phi0,
#                omega_psi=omega_psi, psi0=psi0)

    def _set_topology(self, topology):
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


    def _parameterize(self):
        """Parameterize the potential terms based off the current topology

        Extract the residue identites and atomic indices from the topology in
        order to assign the transferable parameters.
        """

        self._parameterize_debye()
        self._parameterize_rama()
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

#        # pairs included in the local density calculation for each residue. We
#        # exclude -/+1 neighbors on the same chain.
#        self.res_burial_pairs = []
#        for i in range(self.top.n_residues):
#            temp_res_burial_pairs = []
#            resi = self.top.residue(i)
#            if resi.name == "GLY":
#                idx1 = self.top.select("resid {} and name CA".format(i))[0]
#            else:
#                idx1 = self.top.select("resid {} and name CB".format(i))[0]
#            for j in range(self.top.n_residues):
#                if i != j:
#                    resj = self.top.residue(j)
#
#                    if resj.name == "GLY":
#                        idx2 = self.top.select("resid {} and name CA".format(j))[0]
#                    else:
#                        idx2 = self.top.select("resid {} and name CB".format(j))[0]
#
#                    if (resi.chain.index == resj.chain.index):
#                        if abs(resj.index - resi.index) > 1:
#                            temp_res_burial_pairs.append([idx1, idx2])
#                    else:
#                        temp_res_burial_pairs.append([idx1, idx2])
#
#            self.res_burial_pairs.append(np.array(temp_res_burial_pairs))

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
                    if (resj.index - resi.index) >= self._contact_exclude_neighbors:
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


    def _parameterize_debye(self):
        # Only charged residues interact.

        debye_pairs = []
        debye_charges = []
        debye_kij = []
        for i in range(len(self._charged_residues)):
            resi = self.top.residue(self._charged_residues[i][0] - 1)
            for j in range(i + 1, len(self._charged_residues)):
                resj = self.top.residue(self._charged_residues[j][0] - 1)

                if abs(resi.index - resj.index) >= self._debye_exclude_neighbors:
                    qi = self._charged_residues[i][1]
                    qj = self._charged_residues[j][1]
                    idx1 = self.top.select("resid {} and name CB".format(resi.index))[0]
                    idx2 = self.top.select("resid {} and name CB".format(resj.index))[0]
                    debye_pairs.append([idx1, idx2])
                    debye_charges.append([qi, qj])

                    if (qi > 0) and (qj > 0):
                        debye_kij.append(self._debye_kplusplus)
                    elif (qi < 0) and (qj < 0):
                        debye_kij.append(self._debye_kminusminus)
                    elif ((qi > 0) and (qj < 0)) or ((qi < 0) and (qj > 0)):
                        debye_kij.append(self._debye_kplusminus)
        self._debye_pairs = np.array(debye_pairs)
        self._debye_charges = np.array(debye_charges)
        self._debye_kij = np.array(debye_kij)
        self.n_debye_pairs = len(debye_pairs)

    def _parameterize_alpha_helical(self):
        """Parameterize the alpha-helical interaction"""
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

    def _parameterize_rama(self):
        """Parameterize Rama interaction. Proline is treated special"""

        phi_specs = lambda res: [[res.index - 1, "C"], [res.index, "N"], [res.index, "CA"], [res.index, "C"]]
        psi_specs = lambda res: [[res.index, "N"], [res.index, "CA"], [res.index, "C"], [res.index + 1, "N"]]
        select_atm = lambda idx, name: self.top.select("resid {} and name {}".format(idx, name))[0]

        phi_idxs = [] # C_i-1 N_i CA_i C_i
        psi_idxs = [] # N_i CA_i C_i N_i+1
        pro_phi_idxs = []
        pro_psi_idxs = []
        for chain in self.top.chains:
            for res in chain.residues:
                if (res.index in [chain.residue(0).index,
                    chain.residue(chain.n_residues - 1).index]) or (res.name == "GLY"):
                    # We skip terminal residues and glycine
                    continue
                else:
                    res_phi_idxs = [ select_atm(idx, name) for idx, name in phi_specs(res) ]
                    res_psi_idxs = [ select_atm(idx, name) for idx, name in psi_specs(res) ]

                    if res.name == "PRO":
                        # add to proline dihedrals
                        pro_phi_idxs.append(res_phi_idxs)
                        pro_psi_idxs.append(res_psi_idxs)
                    else:
                        # add to regular dihedrals
                        phi_idxs.append(res_phi_idxs)
                        psi_idxs.append(res_psi_idxs)

        self._phi_idxs = np.array(phi_idxs)
        self._psi_idxs = np.array(psi_idxs)
        self._pro_phi_idxs = np.array(pro_phi_idxs)
        self._pro_psi_idxs = np.array(pro_psi_idxs)
        self.n_phi = len(phi_idxs)
        self.n_pro_phi = len(pro_phi_idxs)

#    def _calculate_local_density(self, traj):
#        """Calculate local protein density around each residue"""
#
#        local_density = np.zeros((traj.n_frames, traj.n_residues))
#        direct = self.potential_forms["DIRECT"]
#
#        for i in range(traj.n_residues):
#            # compute local density of residue
#            pairs = self.res_burial_pairs[i]
#            local_density[:,i] = np.sum(direct.theta_I(md.compute_distances(traj, pairs)), axis=1)
#        return local_density

    def _calculate_local_density(self, traj):
        """Calculate local protein density around each residue"""

        local_density = np.zeros((traj.n_frames, traj.n_residues))
        direct = self.potential_forms["DIRECT"]

        is_gly_ca = lambda atom: ((atom.residue.name == "GLY") and (atom.name == "CA"))
        is_other_cb = lambda atom: ((atom.residue.name != "GLY") and (atom.name == "CB"))
        not_neighbors = lambda atom, idxs: not (atom.residue.index in idxs)

        for i in range(traj.n_residues):
            # compute local density of residue
            res = traj.top.residue(i)
            if res.name == "GLY":
                idx1 = traj.top.select("resid {} and name CA".format(res.index))
            else:
                idx1 = traj.top.select("resid {} and name CB".format(res.index))

            # If terminal
            if res.index in self._n_terminal_residues:
                pairs = np.array([ [idx1, atom.index] for atom in traj.top.atoms \
                        if is_gly_ca(atom) and not_neighbors(atom, [res.index, res.index + 1]) ] +\
                        [ [idx1, atom.index] for atom in traj.top.atoms \
                        if is_other_cb(atom) and not_neighbors(atom, [res.index, res.index + 1]) ])
            elif res.index in self._c_terminal_residues:
                pairs = np.array([ [idx1, atom.index] for atom in traj.top.atoms \
                        if is_gly_ca(atom) and not_neighbors(atom, [res.index, res.index - 1]) ] +\
                        [ [idx1, atom.index] for atom in traj.top.atoms \
                        if is_other_cb(atom) and not_neighbors(atom, [res.index, res.index - 1]) ])
            else:
                pairs = np.array([ [idx1, atom.index] for atom in traj.top.atoms \
                        if is_gly_ca(atom) and not_neighbors(atom, [res.index - 1, res.index, res.index + 1]) ] +\
                        [ [ idx1, atom.index] for atom in traj.top.atoms \
                        if is_other_cb(atom) and not_neighbors(atom, [res.index - 1, res.index, res.index + 1]) ])

            local_density[:,i] = np.sum(direct.theta_I(md.compute_distances(traj, pairs)), axis=1)
        return local_density

    def calculate_burial_energy(self, traj, local_density=None, total=True):
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

        if total:
            Vburial = np.zeros(traj.n_frames, float)
        else:
            Vburial = np.zeros((traj.n_frames, self.top.n_residues), float)

        for i in range(self.top.n_residues):
            if total:
                Vburial += burial.V(res_local_density[:,i], self.res_gamma_burial[i])
            else:
                Vburial[:,i] = burial.V(res_local_density[:,i], self.res_gamma_burial[i])

        return Vburial

    def calculate_direct_energy(self, traj, total=True):
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

        if total:
            Vdirect = np.zeros(bb_traj.n_frames, float)
        else:
            Vdirect = np.zeros((bb_traj.n_frames, self.n_pairs), float)

        for i in range(self.n_pairs):
            gamma_direct = self.gamma_direct[self._contact_gamma_idxs[i,0], self._contact_gamma_idxs[i,1]]
            if total:
                Vdirect += direct.V(r[:,i], gamma_direct)
            else:
                Vdirect[:,i] = direct.V(r[:,i], gamma_direct)
        return Vdirect

    def calculate_water_energy(self, traj, local_density=None, total=True, split=False):
        """Calculate the one-body burial potential

        Parameters
        ----------
        traj : mdtraj.Trajectory
            Trajectory to calculate energy over.
        total : opt, bool
            If true (default) return the sum of the burial potentials. If
            false, return the burial energy of each individual residue.
        split : opt, bool
            If False (default) do nothing. If True, return individual energy
            values for water term and protein term.
        """
        water = self.potential_forms["WATER"]

        bb_traj = self.backbone_mapping.map_traj(traj)
        r = md.compute_distances(bb_traj, self._contact_pairs)

        if local_density is None:
            res_local_density = self._calculate_local_density(bb_traj)
        else:
            res_local_density = local_density

        if total:
            Vwater = np.zeros(bb_traj.n_frames, float)
        else:
            Vwater = np.zeros((bb_traj.n_frames, self.n_pairs), float)

        if split:
            split_water = np.zeros((bb_traj.n_frames, self.n_pairs), float)
            split_protein = np.zeros((bb_traj.n_frames, self.n_pairs), float)

        for i in range(self.n_pairs):
            rhoi = res_local_density[:,self._contact_res_idxs[i,0]]
            rhoj = res_local_density[:,self._contact_res_idxs[i,1]]
            gamma_water = self.gamma_water[self._contact_gamma_idxs[i,0], self._contact_gamma_idxs[i,1]]
            gamma_protein = self.gamma_protein[self._contact_gamma_idxs[i,0], self._contact_gamma_idxs[i,1]]
            if not split:
                if total:
                    Vwater += water.V(r[:,i], rhoi, rhoj, gamma_water, gamma_protein)
                else:
                    Vwater[:,i] = water.V(r[:,i], rhoi, rhoj, gamma_water, gamma_protein)
            else:
                split_water[:,i] = -water.lambda_water * gamma_water * water.dVdgamma_water(r[:,i], rhoi, rhoj)
                split_protein[:,i] = -water.lambda_water * gamma_protein * water.dVdgamma_protein(r[:,i], rhoi, rhoj)

        if split:
            return split_water, split_protein
        else:
            return Vwater

    def calculate_debye_energy(self, traj, total=True):
        """Calculate the one-body burial potential

        Parameters
        ----------
        traj : mdtraj.Trajectory
            Trajectory to calculate energy over.
        sum : opt, bool
            If true (default) return the sum of the burial potentials. If
            false, return the burial energy of each individual residue.
        """
        debye = self.potential_forms["DEBYE"]

        bb_traj = self.backbone_mapping.map_traj(traj)

        r = md.compute_distances(bb_traj, self._debye_pairs)

        if total:
            Vdebye = np.zeros(bb_traj.n_frames, float)
        else:
            Vdebye = np.zeros((bb_traj.n_frames, self.n_debye_pairs), float)

        for i in range(self.n_debye_pairs):
            qi = self._debye_charges[i,0]
            qj = self._debye_charges[i,1]
            kij = self._debye_kij[i]
            if total:
                Vdebye += debye.V(r[:,i], qi, qj, kij)
            else:
                Vdebye[:, i] = debye.V(r[:,i], qi, qj, kij)
        return Vdebye

    def calculate_helix_energy(self, traj, local_density=None, total=True):
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

        if total:
            Vhelix = np.zeros(bb_traj.n_frames, float)
        else:
            Vhelix = np.zeros((bb_traj.n_frames, self.n_alpha_helix), float)

        for i in range(self.n_alpha_helix):
            rhoi = res_local_density[:,self._helix_res_idxs[i,0]]
            rhoi_4 = res_local_density[:,self._helix_res_idxs[i,1]]
            fai = self._helix_fai[i]
            fai_4 = self._helix_fai_4[i]
            if total:
                Vhelix += helix.V(r_ON[:,i], r_OH[:,i], rhoi, rhoi_4, fai, fai_4)
            else:
                Vhelix[:, i] = helix.V(r_ON[:,i], r_OH[:,i], rhoi, rhoi_4, fai, fai_4)
        return Vhelix

    def calculate_rama_energy(self, traj, total=True):
        """Calculate the one-body burial potential

        Parameters
        ----------
        traj : mdtraj.Trajectory
            Trajectory to calculate energy over.
        sum : opt, bool
            If true (default) return the sum of the burial potentials. If
            false, return the burial energy of each individual residue.
        """
        rama = self.potential_forms["RAMA"]
        pro_rama = self.potential_forms["RAMA_PROLINE"]

        bb_traj = self.backbone_mapping.map_traj(traj)

        # where does AWSEM define 0.
        phi = md.compute_dihedrals(bb_traj, self._phi_idxs)
        psi = md.compute_dihedrals(bb_traj, self._psi_idxs)

        if total:
            Vrama = np.zeros(bb_traj.n_frames, float)
        else:
            Vrama = np.zeros((bb_traj.n_frames, self.n_phi + self.n_pro_phi), float)

        for i in range(self.n_phi):
            if total:
                Vrama += rama.V(phi[:,i], psi[:,i])
            else:
                Vrama[:, i] = rama.V(phi[:,i], psi[:,i])

        if self.n_pro_phi > 0:
            pro_phi = md.compute_dihedrals(bb_traj, self._pro_phi_idxs)
            pro_psi = md.compute_dihedrals(bb_traj, self._pro_psi_idxs)
            for i in range(self.n_pro_phi):
                if total:
                    Vrama += pro_rama.V(pro_phi[:,i], pro_psi[:,i])
                else:
                    Vrama[:,self.n_phi + i] = pro_rama.V(pro_phi[:,i], pro_psi[:,i])
        return Vrama

    def add_fragment_memory(self, traj, protein_index, frag_index, length, weight):
        #construct list of atoms from the fragment index
        fragment_top = traj.top
        fragment_atom_list = []
        for idx in np.arange(frag_index, frag_index+length):
            for atom in fragment_top.residue(idx).atoms:
                if atom.name in ["CA", "CB"]:
                    fragment_atom_list.append(atom)

        protein_top = self.three_bead_topology #traj files write three beads
        protein_atom_list = []
        for idx in np.arange(protein_index, protein_index+length):
            for atom in protein_top.residue(idx).atoms:
                if atom.name in ["CA", "CB"]:
                    protein_atom_list.append(atom)

        #check the lists, make sure the sequences match
        assert len(protein_atom_list) == len(fragment_atom_list)

        #generate a list of pairs of atoms and indices
        fragment_atom_pairs = []
        distance_pairs = []
        protein_atom_pairs = []
        num_atoms = len(fragment_atom_list)

        for idx in range(num_atoms):
            frag_atm1 = fragment_atom_list[idx]
            prot_atm1 = protein_atom_list[idx]
            for jdx in np.arange(idx+1, num_atoms):
                frag_atm2 = fragment_atom_list[jdx]
                prot_atm2 = protein_atom_list[jdx]
                #only add pairs separated by two residues
                separation = frag_atm1.residue.index-frag_atm2.residue.index
                if np.abs(separation) > 2:
                    fragment_atom_pairs.append([frag_atm1, frag_atm2])
                    protein_atom_pairs.append([prot_atm1, prot_atm2])
                    distance_pairs.append([frag_atm1.index,frag_atm2.index])
                else:
                    pass

        #compute distances
        if np.shape(distance_pairs)[0] == 0:
            raise IOError("No Distance Pairs found!")
        distances = md.compute_distances(traj, distance_pairs, periodic=False)
        distances = distances.transpose()[:,0] * 10.#reform to NX1 array

        #add to the Hamiltonian, each fragment term to a list
        fragment = awsem.AWSEM_POTENTIALS["FRAGMENT"](protein_atom_pairs, distances, weight=weight)

        try:
            self.fragment_potentials.append(fragment)
        except AttributeError:
            self.fragment_potentials = [fragment]

    def calculate_fragment_memory_potential(self, traj, total=True, dgamma=False):
        if not hasattr(self,"fragment_potentials"):
            raise AttributeError("fragment_potentials not initialized")

        energy_list = []
        for idx, potential in enumerate(self.fragment_potentials):
            distances = md.compute_distances(traj, potential.atom_pair_indices, periodic=False) * 10.
            if dgamma:
                energy = potential.dVdgamma(distances)
            else:
                energy = potential.V(distances)

            if idx == 0:
                total_potential = np.copy(energy)
            else:
                total_potential += energy
            energy_list.append(energy)

        if total == False:
            energy_array = np.array(energy_list)
            energy_array *= self.fragment_memory_scale
            return energy_array
        else:
            total_potential *= self.fragment_memory_scale
            return total_potential

    def calculate_energy(self):
        pass

    # TODO:
    # - method to select a subset of interactions to calculate the enegy of.

if __name__ == "__main__":
    pass
