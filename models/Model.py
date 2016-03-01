
import model_builder.models.structure.mapping as mpg
import model_builder.models.potentials as ptl
import model_builder.models.contacts as cts

'''
A Model consists of:
- a structure mapping
- a set of potentials
    - a set of interaction parameters
'''

class Model(object):

    def __init__(self, topology=None, traj=None):
        if (topology is None) and (traj is not None):
            topology = traj.top
        self.structure_mapping = mpg.CalphaMapping(topology)  
        self.potentials = ptl.Potentials()

    def set_reference(self, traj):
        self.ref_traj_aa = traj[0]
        self.ref_traj = self.structure_mapping.map_traj(traj[0])

    def add_sbm_contacts(self):
        residue_contacts = cts.residue_contacts(ref_traj)
        atm_pairs = self.structure_mapping.residue_to_atom_contacts(residue_contacts)

        eps = 1.
        xyz = self.ref_traj.xyz
        self.pairV = []
        for atm1, atm2 in atm_pairs:
            r0 = np.linalg.norm(xyz[atm1.index,:] - xyz[atm2.index,:])
            self.potentials.add_pair(code, atm1, atm2, r0)

    def describe(self):
        pass 

    def Vij(self):
        return [ interaction.Vij for interaction in self.pairV ]

    def map_traj(self, traj):
        self.structure_mapping.map_traj(traj)

