import os

import mdtraj as md

import model_builder as mdb

if __name__ == "__main__":
    # Load all-atom pdb structure.
    name = "SH3"
    traj = md.load(name + ".pdb")

    # Create CA structure-based model.
    model = mdb.models.StructureBasedModel(traj.top, bead_repr="CA")

    # Create the Hamiltonian based off of the reference structure.
    model.set_reference(traj)
    model.add_sbm_potentials()

    # Write the gromacs simulation files.
    writer = mdb.models.output.GromacsFiles(model)
    writer.generate_topology()

    with open("topol.top", "w") as fout:
        fout.write(writer.topfile)

