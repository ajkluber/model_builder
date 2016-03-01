
import model_builder.models.structure.mapping as mapping
import model_builder.models.interactions as interactions
import model_builder.models.contacts as contacts

'''
A Model consists of:
- a structure mapping
- a set of interactions
    - a set of interaction parameters
'''

class Model(object):

    def __init__(self, topology=None, traj=None):
        if (topology is None) and (traj is not None):
            topology = traj.top
        self.structure_mapping = mapping.CalphaMapping(topology)  
        self.potentials = None

    def set_reference(self, traj):
        self.ref_traj_aa = traj[0]
        self.ref_traj = self.structure_mapping.map_traj(traj[0])

    def add_sbm_contacts(self):
        self.pairV = interactions.sbm_contacts(self.structure_mapping, self.ref_traj) 

        residue_contacts = contacts.residue_contacts(ref_traj)
        atm_pairs = self.structure_mapping.residue_to_atom_contacts(residue_contacts)

        xyz = ref_traj[0].xyz
        code = 2
        self.pairV = []
        eps = 1.
        for atm1, atm2 in atm_pairs:
            r0 = np.linalg.norm(xyz[atm1.index,:] - xyz[atm2.index,:])
            self.pairV.append(LJ1210Potential(atm1, atm2, r0, eps))


    def describe(self):
        pass 

    def Vij(self):
        return [ interaction.Vij for interaction in self.pairV ]

    def map_traj(self, traj):
        self.structure_mapping.map_traj(traj)

