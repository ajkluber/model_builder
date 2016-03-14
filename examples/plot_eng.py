import os
import numpy as np 
import subprocess as sb
import matplotlib.pyplot as plt

import mdtraj as md

import model_builder as mdb

def plot_energy(ax, gmxeng, myeng, xlabel, ylabel, title=""):

    maxE = max([max(myeng), max(gmxeng)])
    minE = min([min(myeng), min(gmxeng)])

    # calculate R^2 value

    ax.plot([minE, maxE], [minE, maxE], 'k')
    ax.plot(gmxeng, myeng, 'r.')
    ax.set_ylim(minE, maxE)
    ax.set_xlim(minE, maxE)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_title(title, fontsize=18)


if __name__ == "__main__":
    # Load all-atom pdb structure.
    name = "SH3"
    reftraj = md.load(name + ".pdb")

    # Create CA structure-based model.
    model = mdb.models.StructureBasedModel(reftraj.top, bead_repr="CA")

    # Create the Hamiltonian based off of the reference structure.
    model.set_reference(reftraj)
    model.add_sbm_potentials()

    traj = md.load("traj.xtc", top=model.mapping.top)

    if not os.path.exists("Eterms.xvg"):
        sb.call("""g_energy_sbm -f ener.edr -s topol.tpr -xvg none -o Eterms.xvg << HERE
Bond
Angle
Proper-Dih.
LJ-14
LJ-(SR)
Potential
Total-Energy
HERE""", shell=True)

    Eterms = np.loadtxt("Eterms.xvg")
    
    # calculate energy terms using model_builder
    Ebond = model.Hamiltonian.calc_bond_energy(traj)
    Eangle = model.Hamiltonian.calc_angle_energy(traj)
    Edih = model.Hamiltonian.calc_dihedral_energy(traj)
    Epair = model.Hamiltonian.calc_pair_energy(traj)

    # plot comparison
    fig, axes = plt.subplots(2, 2, figsize=(14,10))

    plot_energy(axes[0,0], Eterms[:,1], Ebond, "gmx bond", "mdb bond")
    plot_energy(axes[0,1], Eterms[:,2], Eangle, "gmx angle", "mdb angle")
    plot_energy(axes[1,0], Eterms[:,3], Edih, "gmx dihedral", "mdb dihedral")
    plot_energy(axes[1,1], Eterms[:,4], Epair, "gmx pair", "mdb pair")

    fig.suptitle("Energy terms", fontsize=20)

    if not os.path.exists("plots"):
        os.mkdir("plots")
    os.chdir("plots")
    fig.savefig("compare_energy_terms.pdf", bbox_inches="tight")
    fig.savefig("compare_energy_terms.png", bbox_inches="tight")
    os.chdir("..")
