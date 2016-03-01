import model_builder.models.structure.contacts as contacts

import model_builder.models.interactions.pair_potentials as pair_potentials

class PairPotential(object):
    
    def __init__(self, atmi, atmj):
        self.atmi = atmi
        self.atmj = atmj

    def __repr__(self):
        """interaction description"""
        return "%s:%12s%12s" % (self.prefix_label, self.atmi, self.atmj)

class LJ1210Potential(PairPotential):
    
    def __init__(self, atmi, atmj, r0, eps):
        PairPotential.__init__(self, atmi, atmj)
        #self.code = 2
        self.prefix_label = "LJ1210"
        self.r0 = r0
        self.eps = eps

    def V(self, r):
        return self.eps*self.dim_V(r)

    def dim_V(self, r):
        x = self.r0/r
        V = 5.*(x**12) - 6.*(x**10)
        return V

    def dVdr(self, r):
        return self.eps*self.dim_dVdr(r)

    def dim_dVdr(self, r):
        x = self.r0/r
        V = (-60./self.r0)*((x**13) - (x**11))
        return V
    
def sbm_contacts(mapping, ref_traj):
    residue_contacts = contacts.residue_contacts(ref_traj)
    atm_pairs = mapping.residue_to_atom_contacts(residue_contacts)

    xyz = ref_traj[0].xyz 
    code = 2
    pairV = []
    eps = 1.
    for atm1, atm2 in atm_pairs:
        r0 = np.linalg.norm(xyz[atm1.index,:] - xyz[atm2.index,:])
        pairV.append(LJ1210Potential(atm1, atm2, r0, eps))
    return pairV

if __name__ == "__main__":
    pass
