import os
import numpy as np 
import subprocess as sb
import matplotlib.pyplot as plt

import mdtraj as md

import model_builder as mdb

import run_butane

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
    save = True

    xyz, top, model = run_butane.butane_toy_model()

    os.chdir("sim_data")
    # Get simulation trajectory
    traj = md.load("traj.xtc", top=top)

    # Get energy from simulation
    Ebond_gmx = np.loadtxt("Ebond_gmx.xvg", usecols=(,1))
    Eangle_gmx = np.loadtxt("Eangle_gmx.xvg", usecols=(,1))
    Edih_gmx = np.loadtxt("Edih_gmx.xvg", usecols=(,1))
    os.chdir("..")
    
    # calculate energy terms using model_builder
    Ebond_mdb = model.Hamiltonian.calc_bond_energy(traj)
    Eangle_mdb = model.Hamiltonian.calc_angle_energy(traj)
    Edih_mdb = model.Hamiltonian.calc_dihedral_energy(traj)

    # plot comparison of energy terms
    fig, axes = plt.subplots(2, 2, figsize=(14,10))

    plot_energy(axes[0,0], Ebond_gmx, Ebond_mdb, "gmx bond", "mdb bond")
    plot_energy(axes[0,1], Eangle_gmx, Eangle_mdb, "gmx angle", "mdb angle")
    plot_energy(axes[1,0], Edih_gmx, Edih_mdb, "gmx dihedral", "mdb dihedral")

    fig.suptitle("Energy terms", fontsize=20)

    if save:
        if not os.path.exists("plots"):
            os.mkdir("plots")
        os.chdir("plots")
        fig.savefig("compare_energy_terms.pdf", bbox_inches="tight")
        fig.savefig("compare_energy_terms.png", bbox_inches="tight")
        os.chdir("..")

    plt.show()
