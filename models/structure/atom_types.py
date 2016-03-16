
class CoarseGrainAtom(object):
    def __init__(self, index, name, residx, resname, radius, mass, charge):
        self.index = index
        self.name = name
        self.residx = residx
        self.resname = resname
        self.radius = radius
        self.mass = mass
        self.charge = charge
        self.ptype = "A"    # Atom (not virtual site)

        # Lennard-Jones parameters
        self.c6 = 0
        self.c12 = (self.radius)**12

    def describe(self):
        return "<model_builder.CoarseGrainAtom {} {} {} {}>".format(
                self.name, self.radius, self.mass, self.charge)

