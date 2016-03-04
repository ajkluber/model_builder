


class CoarseGrainAtom(object):
    def __init__(self, name, radius, mass, charge):
        self.name = name
        self.radius = radius
        self.mass = mass
        self.charge = charge

    def describe(self):
        return "<model_builder.CoarseGrainAtom {} {} {} {}>".format(
                self.name, self.radius, self.mass, self.charge)

