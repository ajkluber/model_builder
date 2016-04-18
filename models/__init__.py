'''
Classes
-------
Model

StructureBasedModel

'''

import numpy as np

import output

import structure.mappings as mpg
import structure.contacts as cts
import potentials as ptl


class Model(object):
    """Model class """
    def __init__(self, topology, bead_repr="CA"):
        self.mapping = mpg.assign_mapping(bead_repr, topology) 
        self.Hamiltonian = ptl.Hamiltonian()

    def describe(self):
        #TODO: What is the best description?
        pass 
    
    def map_traj(self, traj):
        return self.mapping.map_traj(traj)

    def save_starting_conf(self, saveas="conf.gro"):
        self.starting_traj.save(saveas)

    def set_starting_conf(self, traj):
        self.starting_traj = traj
        
    def add_pairs(self, pairs):
        self.mapping._add_pairs(pairs)
    
    def assign_fitted_epsilons(self, params_to_fit_indices):
        self.params_to_fit_indices = params_to_fit_indices
        self.fitted_epsilons = []
        for i in params_to_fit_indices:
            self.fitted_epsilons.append(self.Hamiltonian._epsilons[i])   
    
    def output_epsilons(self):
        params = self.Hamiltonian._epsilons
        for idx, i in enumerate(self.params_to_fit_indices):
            params[i] = self.fitted_epsilons[idx]
        
        f = open("params", "w")
        f.write("# Fitted Parameters\n")
        for param in params:
            f.write("%f\n"%param)
        f.close()
            
        
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
        self.Hamiltonian = ptl.StructureBasedHamiltonian()
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
    
    def assign_disulfides(self, disulfides):
        self.mapping.add_disulfides(disulfides)
        
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
    
    
    
    
    
        
