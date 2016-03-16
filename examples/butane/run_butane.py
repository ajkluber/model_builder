import os
import numpy as np
import subprocess as sb
import matplotlib.pyplot as plt

import mdtraj as md
from mdtraj.core.element import get_by_symbol

import model_builder as mdb

import simulation.mdp

#def generate_dihedral_traj():
#    n_points = 1000
#    phi = np.linspace(0, 2.*np.pi, n_points)
#
#    xyz = np.zeros((n_points, 4, 3), float) 
#    xyz[:,0,:] = np.array([1, -1, 0])
#    xyz[:,1,:] = np.array([0, -1, 0])
#    xyz[:,2,:] = np.array([0, 0, 0])
#    xyz[:,3,0] = np.cos(phi)
#    xyz[:,3,2] = np.sin(phi)
#    
#    traj = md.Trajectory(xyz, top)
#
#    dih_idxs = np.array([[0, 1, 2, 3]])
#
#    dih = md.compute_dihedrals(traj, dih_idxs)
#
#    plt.plot(phi, dih)
#    plt.plot(phi, np.pi - dih)
#    plt.show()

def butane_toy_model():

    # Starting configuration is close to a cis conformation
    phi0 = np.pi
    phi_ini = 2*np.pi*0.8
    xyz = np.zeros((1, 4, 3), float)
    xyz[:,0,:] = np.array([1, -1, 0])
    xyz[:,1,:] = np.array([0, -1, 0])
    xyz[:,2,:] = np.array([0, 0, 0])
    xyz[:,3,0] = np.cos(phi_ini)
    xyz[:,3,2] = np.sin(phi_ini)

    # Coordinates must be positive in gromacs
    shift = np.array([5,5,5])
    xyz += shift

    # Create a mdtraj Topology for our butane molecule.
    newtop = md.Topology()
    chain = newtop.add_chain()
    for i in range(4):
        res = newtop.add_residue("GLY", chain, i)
        new_ca = newtop.add_atom('CA', get_by_symbol('C'), res, serial=i)
        #if i >= 1:
        #    prev_ca = chain.atom(i - 1)
        #    newtop.add_bond(prev_ca, new_ca)
    model = mdb.models.Model(newtop, bead_repr="CA")

    model.mapping.add_atomtypes()

    top = model.mapping.top

    # Add bond interactions
    model.Hamiltonian._add_bond("HARMONIC_BOND", top.atom(0), top.atom(1), 1, 1)
    model.Hamiltonian._add_bond("HARMONIC_BOND", top.atom(1), top.atom(2), 1, 1)
    model.Hamiltonian._add_bond("HARMONIC_BOND", top.atom(2), top.atom(3), 1, 1)

    # Add angle interactions
    model.Hamiltonian._add_angle("HARMONIC_ANGLE", top.atom(0), top.atom(1),
                                    top.atom(2), 2*(np.pi/180)**2, 90)

    model.Hamiltonian._add_angle("HARMONIC_ANGLE", top.atom(1), top.atom(2),
                                    top.atom(3), 2*(np.pi/180)**2, 90)

    # Add dihedral interaction
    model.Hamiltonian._add_dihedral("COSINE_DIHEDRAL", top.atom(0), top.atom(1),
                                    top.atom(2), top.atom(3), 1, phi0, 1)

    #model.Hamiltonian._add_dihedral("COSINE_DIHEDRAL", top.atom(0), top.atom(1),
    #                                top.atom(2), top.atom(3), 0.5, phi0, 3)

    return xyz, top, model

if __name__ == "__main__":

    # Create
    xyz, top, model = butane_toy_model()

    traj = md.Trajectory(xyz, top)

    if os.path.exists("sim_data"):
        os.rmdir("sim_data")
    os.mkdir("sim_data")
    os.chdir("sim_data")

    # We're going to run a constant energy (newtonian dynamics) simulation to
    # see the oscillations in energy.
    with open("run.mdp", "w") as fout:
        fout.write(simulation.mdp.constant_energy("1000000"))

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
    sb.call("""g_energy_sbm -f ener.edr -s topol.tpr -xvg none -o Epot_gmx.xvg << HERE
Potential
HERE""", shell=True)

    os.chdir("..")
