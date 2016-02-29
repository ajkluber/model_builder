
import model_builder.structure.mapping 
import model_builder.interactions

'''
A Model consists of:
- a structure mapping
- a set of interactions
    - a set of interaction parameters
'''

class Model(object):

    def __init__(self, topology):
        self.structure_mapping = model_builder.structure.mapping.CalphaMapping(topology)  

    def set_reference(self, traj):
        self.aa_ref_xyz = traj[0].xyz 
        self.ref_xyz = self.structure_mapping.map_traj(traj[0]).xyz

    def add_contacts(self):
        self.structure_mapping

    def map_traj(self, traj):
        self.structure_mapping.map_traj(traj)

