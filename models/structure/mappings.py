import mdtraj as md


class StructureMapping(object):

    def __init__(self):
        pass
        

class CalphaMapping(object):

    def __init__(self, traj):
        top = traj.topology
        self.ca_idxs = np.array([ atm.index for atm in top.atoms if atm.name == "CA" ])

    def map(self, traj):
        return traj.xyz[:,self.ca_idxs,:]

class CalphaCbetaMapping(object):

    def __init__(self, traj):
        top = traj.topology
        self.ca_idxs = np.array([ atm.index for atm in top.atoms if atm.name == "CA" ])

    def map(self, traj):
        return traj.xyz[:,self.ca_idxs,:]
