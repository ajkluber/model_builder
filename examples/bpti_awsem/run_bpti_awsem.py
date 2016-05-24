import os
import subprocess as sb

import mdtraj as md

import model_builder as mdb
import simulation

if __name__ == "__main__":
    # The default parameters are in a subdirectory where you downloaded the
    # AWSEM source code. The LAMMPS executable.
    awsem_param_path = "/home/alex/packages/awsemmd/parameters"   
    lammps_exe = "/home/alex/packages/lammps-9Oct12/src/lmp_serial"

    # Creating the AWSEM model can take a while (~30 secs for 58 residues).
    name = "5pti"
    disulfides = [[5, 55],[14, 38],[30, 51]]

    traj = md.load("{}.pdb".format(name))
    model = mdb.models.AwsemModel(traj.top)
    model.source_parameters(awsem_param_path)
    model.set_starting_conf(model.map_traj(traj))
    model.mapping.add_disulfides(disulfides)

    # The two ambiguous energy scales in AWSEM are the temperature and the
    # strength of the "fragment memory" structural bias. Set frag_strength =
    # None if you don't want to use the fragment memory interaction.
    T = 300.
    frag_strength = 0.05
    
    n_frames = int(1E5)
    n_steps_out = 100
    n_steps = n_frames*n_steps_out

    # Structural bias towards a reference structure. Called using 'single
    # memory' AWSEM.
    if not os.path.exists("fraglib"):
        os.mkdir("fraglib")
    model.starting_traj.save("fraglib/{}.gro".format(name))
    start_idx = 1
    frag_mem_string = "[Memories]\n"
    for i in range(native_traj.n_chains):
        n_res = native_traj.top.chain(i).n_residues
        frag_mem_string += "../fraglib/{}.gro {} {} {} 1\n".format(name, start_idx, start_idx, n_res - 1)
        start_idx += n_res

    if not os.path.exists("nvt_data"):
        os.mkdir("nvt_data")

    # Simulate at constant temperature using Langevin dynamics.
    simulation.lammps.awsem.prep_constant_temp(
            model, traj, name, T, n_steps, n_steps_out, frag_strength,
            frag_mem_string, awsem_param_path=awsem_param_path)

    # Run 
    sb.call("{} < {}.in".format(lammps_exe, name), shell=True)

    os.chdir("nvt_data")
