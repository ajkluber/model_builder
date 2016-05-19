

def theta_I(r, nu=50., r_min=0.45, r_max=0.65):
    return 0.25*(1. + np.tanh(nu*(r - r_min)))*(1. + np.tanh(nu*(r_max - r)))


class DirectContact(object):

    def __init__(self, gamma_direct, nu=0.5, r_min=0.45, r_max=0.65):
        self.gamma_direct = gamma_direct
        self.nu = nu
        self.r_min = r_min
        self.r_max = r_max

    def V(self, r): 
        return self.gamma_direct*self.theta_I(r)

    def dVdgamma_direct(self, r):
        return self.theta_I(r)

    def theta_I(self, r):
        return 0.25*(1. + np.tanh(self.nu*(r - self.r_min)))*(1. + np.tanh(self.nu*(self.r_max - r)))

class WaterMediatedContact(object):

    def __init__(self, gamma_water, gamma_protein, nu=0.5, nu_sigma=0.7, r_min=0.65, r_max=0.95):
        self.gamma_water = gamma_water
        self.gamma_protein = gamma_protein
        self.nu = nu
        self.nu_sigma = nu_sigma
        self.r_min = r_min
        self.r_max = r_max

    def V(self, r, rhoi, rhoj): 
        return self.gamma_water*self.dVdgamma_water(r, rhoi, rhoj) +\
            self.gamma_protein*self.dVdgamma_protein(r, rhoi, rhoj)

    def dVdgamma_water(self, r, rhoi, rhoj):
        return self.theta_II(r)*self.sigma_water(rhoi, rhoj)

    def dVdgamma_protein(self, r, rhoi, rhoj):
        return self.theta_II(r)*(1. - self.sigma_water(rhoi, rhoj))

    def theta_II(self, r):
        return 0.25*(1. + np.tanh(self.nu*(r - self.r_min)))*(1. + np.tanh(self.nu*(self.r_max - r)))

    def sigma_water(self, rhoi, rhoj):
        return 0.25*(1. - np.tanh(self.nu_sigma*(rhoi - self.rho_0)))*(1. - np.tanh(self.nu_sigma*(rhoj - self.rho_0)))
    
class Burial(object):
    """Technically a one-body potential"""
    def __init__(self, atmi):
        self.atmi = atmi

    def V(self, r, rhoi):
        pass

AWSEM_POTENTIALS = {"DIRECT":DirectContact,
                "WATER":WaterMediatedContact,
                "BURIAL":Burial}
