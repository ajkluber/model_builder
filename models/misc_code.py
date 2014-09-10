
"""
A snippet excised from SmogCalpha to call shadow.jar externally.

Decided to make the contacts file mandatory input instead.

Also a snippet that gets the Hydrogen bonds

"""

import os
import numpy as np
import mdtraj as md

def shadow_contacts(self):
    """
     Call SMOG Shadow jar code to determine the shadow contacts. If 
        the reference matrix Qref_cryst.dat doesn't exist then create 
        and dive into a subdirectory to run shadow map. Then save 
        Qref_cryst.dat in the parent directory.
    """

    

    subdir = self.subdir

    cwd = os.getcwd()
    print "  Calculating native contacts for: ",subdir
    if os.path.exists(cwd+"/"+subdir+"/contacts.dat"):
        print "  Native contact map "+subdir+"/Qref_cryst.dat exists."
        print "  Skipping shadow map calculation. Loading native contact map..."
        self.contacts = np.loadtxt(cwd+"/"+subdir+"/contacts.dat",dtype=int)
    else:
        print "  Native contact map "+subdir+"/Qref_cryst.dat does not exist."
        print "  Doing shadow map calculation..."
        print "  *** NOTE: module load jdk/1.7.0.21 required for shadow map ***"
        #if not os.path.exists("Qref_shadow"):
        os.chdir(cwd+"/"+subdir+"/Qref_shadow")
        cmd0 = 'cp /projects/cecilia/SCM.1.31.jar .'
        sb.call(cmd0,shell=True,stdout=open("contacts.out","w"),stderr=open("contacts.err","w"))
        cmd1 = 'echo -e "6\\n6\\n" | pdb2gmx -f clean.pdb -o %s.gro -p %s.top -ignh' % (subdir,subdir)
        sb.call(cmd1,shell=True,stdout=open("convert.out","w"),stderr=open("convert.err","w"))
        cmd2 = 'java -jar SCM.1.31.jar -g %s.gro -t %s.top -o %s.contacts -m shadow -c 4.5 --coarse CA' % (subdir,subdir,subdir)
        sb.call(cmd2,shell=True,stdout=open("contacts.out","w"),stderr=open("contacts.err","w"))
        self.contacts = np.loadtxt(subdir+".contacts",dtype=int,usecols=(1,3))
        print "  Native contact map calculated with shadow map. Saving Qref_cryst.dat..."
        np.savetxt(cwd+"/"+subdir+"/contacts.dat",self.contacts,delimiter=" ",fmt="%4d")

    Qref = np.zeros((self.n_residues,self.n_residues))
    for pair in self.contacts:
        Qref[pair[0]-1,pair[1]-1] = 1

    np.savetxt(cwd+"/"+subdir+"/Qref_cryst.dat",Qref,delimiter=" ",fmt="%1d")
    os.chdir(cwd)
    print "  Length = %d  Number of contacts = %d  Nc/L=%.4f" % (len(Qref),sum(sum(Qref)),float(sum(sum(Qref)))/float(len(Qref)))
    self.Qref = Qref
    self.n_contacts = len(self.contacts)

def extract_backbone_Hbonds(self,System):
    """ Use Gromacs and MDTraj to get backbone hydrogen bonds """

    name = System.subdir

    cwd = os.getcwd()
    
    if not os.path.exists(System.path+"/"+name+"/hbonds"):
        os.makedirs(System.path+"/"+name+"/hbonds")
    if not os.path.exists(System.path+"/"+name+"/hbonds/hbonds.dat"):
        ## Determine backbone H-bonds with Gromacs and MDTraj
        shutil.copy(name+".pdb",System.path+"/"+name+"/hbonds/")
        os.chdir(System.path+"/"+name+"/hbonds")
        prep_script = "#!/bin/bash\n"
        prep_script += "echo -e '1\\n6\\n' | pdb2gmx -f %s.pdb -o %s.gro\n" % (name,name)
        prep_script += "trjconv -f %s.gro -o %s.xtc\n" % (name,name)
        prep_script += "editconf -f %s.gro -o %s_withH.pdb\n" % (name,name)
        open("prep_pdb.sh","w").write(prep_script)
        sb.call("bash prep_pdb.sh", shell=True,stdout=open("hbonds.out","w"),stderr=open("hbonds.err","w"))

        traj = md.load("%s.xtc" % name, top="%s_withH.pdb" % name)
        M = (md.kabsch_sander(traj))[0]
        Hbonds = np.array(M.todense())
        np.savetxt("hbonds.dat",Hbonds)
        os.chdir("..")
    else:
        Hbonds = np.loadtxt(System.path+"/"+name+"/hbonds/hbonds.dat")
    os.chdir(cwd)

    self.Hbonds = Hbonds
