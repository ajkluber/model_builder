import os
import numpy as np
import subprocess as sb

import mdtraj as md

import model_builder as mdb
import model_builder.models.potentials as ptl
import simulation.mdp
import plot_eng

class HeavyStructureBasedModel(mdb.models.StructureBasedModel):

    def __init__(self, topology, mass, bead_repr=None):
        """Structure-based Model (SBM)

        Parameters
        ----------
        topology : mdtraj.Topology object
            An mdtraj Topology object that describes the molecular topology.

        bead_repr : str [CA, CACB]
            A code specifying the desired coarse-grain mapping. The all-atom 
        to coarse-grain mapping.

        """

        mdb.models.Model.__init__(self, topology, bead_repr=bead_repr)
        self.Hamiltonian = ptl.StructureBasedHamiltonian()
        self.mapping.add_atoms(mass=mass)

        self._default_parameters = {"kb":5., # kJ/(mol nm^2)
                                    "ka":2.*((np.pi/180.)**2),  # kJ/(mol deg^2)
                                    "kd":1.,     # kJ/mol
                                    "eps":1}     # kJ/(mol nm)


def sh3_ca_sbm():
    # Load all-atom pdb structure.
    name = "SH3"
    traj = md.load(name + ".pdb")

    # Create CA structure-based model.
    model = HeavyStructureBasedModel(traj.top, 100, bead_repr="CA")
    #model = mdb.models.StructureBasedModel(traj.top, bead_repr="CA")

    # Create the Hamiltonian based off of the reference structure.
    model.set_reference(traj)
    model.add_sbm_potentials()

    return model

def save_top_file(model, traj):

    # Write the Hamiltonian Gromacs input file: topol.top
    writer = mdb.models.output.GromacsFiles(model)
    writer.generate_topology()

    with open("topol.top", "w") as fout:
        fout.write(writer.topfile)

    with open("index.ndx", "w") as fout:
        fout.write(writer.index_ndx)

    np.savetxt("table.xvg", writer.tablep, fmt="%16.15e")
    np.savetxt("tablep.xvg", writer.tablep, fmt="%16.15e")
    
    # Save the starting configuration, but we have to fix the unitcell
    # information.
    traj[0].save("conf.gro")
    with open("conf.gro", "r") as fin:
        temp = reduce(lambda x,y: x+y, fin.readlines()[:-1])
        temp += "{:>10f}{:>10f}{:>10f}\n".format(10,10,10)

    with open("conf.gro", "w") as fout:
        fout.write(temp)

    # Save a pdb with conect directives for ease of visualization.
    traj[0].save("ca.pdb")
    mdb.models.structure.viz_bonds.write_bonds_conect(model.mapping.top)
    with open("ca.pdb", "r") as fin:
        temp = reduce(lambda x,y: x+y, fin.readlines()[:-3])
        temp += open("conect.pdb","r").read() + "END\n"
    with open("ca.pdb", "w") as fout:
        fout.write(temp)


def run_and_analyze():
    # Run the simulation
    sb.call("grompp_sbm -f run.mdp -c conf.gro -p topol.top -o topol.tpr", shell=True)
    sb.call("mdrun_sbm -s topol.tpr -table table.xvg -tablep tablep.xvg", shell=True)

    # Calculate each energy term
    sb.call("""g_energy_sbm -f ener.edr -s topol.tpr -xvg none -o Ebond_gmx.xvg << HERE
Bond
HERE""", shell=True)
    sb.call("""g_energy_sbm -f ener.edr -s topol.tpr -xvg none -o Eangle_gmx.xvg << HERE
Angle
HERE""", shell=True)
    sb.call("""g_energy_sbm -f ener.edr -s topol.tpr -xvg none -o Edih_gmx.xvg << HERE
Proper-Dih.
HERE""", shell=True)
    sb.call("""g_energy_sbm -f ener.edr -s topol.tpr -xvg none -o Epair_gmx.xvg << HERE
LJ-14
HERE""", shell=True)
    sb.call("""g_energy_sbm -f ener.edr -s topol.tpr -xvg none -o Epot_gmx.xvg << HERE
Potential
HERE""", shell=True)

if __name__ == "__main__":
    model = sh3_ca_sbm()
    startingtraj = model.ref_traj

    ##################################
    # Run a constant Energy simlation
    ##################################
    if not os.path.exists("nve_data"):
        os.mkdir("nve_data")
    os.chdir("nve_data")

    # Save mdp file for NVE simulation (newtonian dynamics). We want to see the
    # oscillations in energy.
    with open("run.mdp", "w") as fout:
        fout.write(simulation.mdp.constant_energy("1000000"))

    # Save remaining simulation files
    save_top_file(model, startingtraj)
    run_and_analyze()
    plot_eng.plot_energy_terms(model, display=True)
    os.chdir("..")


    #######################################
    # Run a constant temperature simulation
    #######################################
    if not os.path.exists("nvt_data"):
        os.mkdir("nvt_data")
    os.chdir("nvt_data")

    # Save mdp file for NVT simulation (constant temperature).
    with open("run.mdp", "w") as fout:
        fout.write(simulation.mdp.constant_temperature(10, "1000000"))

    # Save remaining simulation files
    save_top_file(model, startingtraj)
    run_and_analyze()
    plot_eng.plot_energy_terms(model, display=True)

    os.chdir("..")
