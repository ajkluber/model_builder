import os
import numpy as np 
import subprocess as sb
import matplotlib.pyplot as plt

import mdtraj as md

import model_builder as mdb

def plot_energy(ax, data1, data2, xlabel, ylabel, title="", ls='', withxy=False):

    if withxy:
        maxdata = max([max(data2), max(data1)])
        mindata = min([min(data2), min(data1)])
        ax.plot([mindata, maxdata], [mindata, maxdata], 'k')
        ax.set_ylim(mindata, maxdata)
        ax.set_xlim(mindata, maxdata)

        # calculate R^2 value
        r2 = np.mean((data1 - data1.mean())*(data2 - data2.mean()))/(np.std(data1)*np.std(data2))
        ax.annotate("$R^2 = {:.5f}$".format(r2), xy=(0,0), xytext=(0.1,0.8), xycoords="axes fraction",
                    textcoords="axes fraction",fontsize=18)

    if ls == '':
        ax.plot(data1, data2)
    else:
        ax.plot(data1, data2, ls)

    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_title(title, fontsize=18)

def plot_energy_terms(model, save=True, display=False):
    # Get simulation trajectory
    traj = md.load("traj.xtc", top=model.mapping.top)

    # Get energy from simulation
    t = np.loadtxt("Ebond_gmx.xvg", usecols=(0,))
    Ebond_gmx = np.loadtxt("Ebond_gmx.xvg", usecols=(1,))
    Eangle_gmx = np.loadtxt("Eangle_gmx.xvg", usecols=(1,))
    Edih_gmx = np.loadtxt("Edih_gmx.xvg", usecols=(1,))
    Epot_gmx = np.loadtxt("Epot_gmx.xvg", usecols=(1,))
    Epair_gmx = np.loadtxt("Epair_gmx.xvg", usecols=(1,))
    
    # calculate energy terms using model_builder
    Ebond_mdb = model.Hamiltonian.calc_bond_energy(traj)
    Eangle_mdb = model.Hamiltonian.calc_angle_energy(traj)
    Edih_mdb = model.Hamiltonian.calc_dihedral_energy(traj)
    Epair_mdb = model.Hamiltonian.calc_pair_energy(traj)

    # Plot correlation of model_builder and gromacs
    fig, axes = plt.subplots(2, 2, figsize=(14,10))

    plot_energy(axes[0,0], Ebond_gmx, Ebond_mdb, "GMX $E_{bond}$", "MDB $E_{bond}$", ls='r.', withxy=True)
    plot_energy(axes[0,1], Eangle_gmx, Eangle_mdb, "GMX $E_{angle}$", "MDB $E_{angle}$", ls='r.', withxy=True)
    plot_energy(axes[1,0], Edih_gmx, Edih_mdb, "GMX $E_{dih}$", "MDB $E_{dih}$", ls='r.', withxy=True)
    plot_energy(axes[1,1], Epair_gmx, Epair_mdb, "GMX $E_{pair}$", "MDB $E_{pair}$", ls='r.', withxy=True)

    # Plot model_builder and gromacs timeseries
    fig2, axes2 = plt.subplots(2, 2, figsize=(14,10))

    plot_energy(axes2[0,0], t, Ebond_gmx, "time", "$E_{bond}$")
    plot_energy(axes2[0,0], t, Ebond_mdb, "time", "$E_{bond}$", ls='--')

    plot_energy(axes2[0,1], t, Eangle_gmx, "time", "$E_{angle}$")
    plot_energy(axes2[0,1], t, Eangle_mdb, "time", "$E_{angle}$", ls='--')

    plot_energy(axes2[1,0], t, Edih_gmx, "time", "$E_{dih}$")
    plot_energy(axes2[1,0], t, Edih_mdb, "time", "$E_{dih}$", ls='--')

    plot_energy(axes2[1,1], t, Epair_gmx, "time", "$E_{pair}$")
    plot_energy(axes2[1,1], t, Epair_mdb, "time", "$E_{pair}$", ls='--')

    if save:
        if not os.path.exists("plots"):
            os.mkdir("plots")
        os.chdir("plots")
        fig.savefig("Eterms_corr.pdf", bbox_inches="tight")
        fig.savefig("Eterms_corr.png", bbox_inches="tight")
        fig2.savefig("Eterms_vs_time.pdf", bbox_inches="tight")
        fig2.savefig("Eterms_vs_time.png", bbox_inches="tight")
        os.chdir("..")

    if display:
        plt.show()

if __name__ == "__main__":
    pass
