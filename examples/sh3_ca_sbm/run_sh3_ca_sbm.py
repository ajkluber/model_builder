import os
import numpy as np
import subprocess as sb

import mdtraj as md

import model_builder as mdb
import simulation
import plot_eng

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
    model, fitopts = mdb.inputs.load_model("SH3.ini")
    writer = mdb.models.output.GromacsFiles(model)

    ##################################
    # Run a constant Energy simlation
    ##################################
    if not os.path.exists("nve_data"):
        os.mkdir("nve_data")
    os.chdir("nve_data")

    # Save mdp file for NVE simulation (newtonian dynamics). We want to see the
    # oscillations in energy.
    with open("run.mdp", "w") as fout:
        fout.write(simulation.gromacs.mdp.constant_energy("1000000"))

    # Save remaining simulation files
    writer.write_simulation_files()
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
        fout.write(simulation.gromacs.mdp.constant_temperature(10, "1000000"))

    # Save remaining simulation files
    writer.write_simulation_files()
    run_and_analyze()
    plot_eng.plot_energy_terms(model, display=True)

    os.chdir("..")
