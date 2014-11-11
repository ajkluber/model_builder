Setting Up
==========

The info below is how to setup and use this package on the Davinci cluster at
Rice University (however majority of info is system independent).


Dependecies
-----------

Here is the fast version of the dependencies. Below is step-by-step 
instructions.

1. python2.7
2. numpy
3. scipy
4. matplotlib
5. ipython
6. pandas
7. Cython
8. nose
9. mdtraj    <-- Package for analyzing MD trajectories
10. mpi4py
11. modeller <-- Package for creating PDB models. Adding loops/mutations.



Installation overview
---------------------

Davinci has python2.6 by default, which comes with several important packages
(namely numpy and matplotlib). We would like a python environment with the most
up-to-date and reliable packages (with the proper configuration).

I created ~/packages to put the tarballs (.tar.gz files) and perform the build
and installation process for each package. We also need to make the ~/.local
directory for our local installation. 

The trickiest and most important package is numpy. Because many other packages
depend on numpy internally, messing up the numpy installation can lead to
strange, cryptic errors with other packages down the road. Also because numpy
is a common dependency it should be the first package installed.

Almost every python package can be installed in one or two painless lines.
Sometimes "python setup.py install --user" will be enough, but numpy is
different because it should be compiled with Intel MKL optimized routines for
best performance. On Davinci, the Intel MKL libaries are found in
/opt/apps/intel/2013.1.039/mkl (NOTE: This location changes with time!!), but
the last part of the path (with the numbers) may change when they update the
server.

After installing python2.7, numpy, and scipy the order doesn't matter, but I
retain the order I did everything for accuracy. All packages except two are
standard python packages. The two third-party packages are mdtraj, developed in
Vijay Pande's group at Stanford and MODELLER developed in Andre Sali's lab at
University of California San Francisco. MODELLER installation requires that you
request a registration key. They can be accessed at the following links:

`mdtraj <http://mdtraj.org>`_
`MODELLER  <http://salilab.org/modeller>`_

The rest of the python packages can be found online in the Python Package Index
PyPI.

List and order of packages installed in this HOW-TO:

1. python2.7
2. numpy
3. scipy
4. matplotlib
5. ipython
6. pandas
7. Cython
8. nose
9. mdtraj    <-- Package for analyzing MD trajectories
10. mpi4py
11. modeller <-- Package for creating PDB models. Adding loops/mutations.


Instructions
------------
1. python2.7 
^^^^^^^^^^^^
Directions sourced from `here <http://isezen.com/2011/09/02/how-to-install-locally-python-on-linux-home-directory>`_

    
Download python2.7 from `here <https://www.python.org/download/releases/2.7.6>`_ 
and put it in the ~/packages directory, then untar and enter::

    cd ~/packages
    tar xzvf Python-2.7.6.tgz
    cd Python-2.7.6

Configure for local installation::

    ./configure
    make altinstall prefix=~/.local exec-prefix=~/.local 

Create a link for convenience::

    cd ~/.local/bin
    ln -s python2.7 python


For convenience define the following alias in your ~/.bashrc (Remember to
source ~/.bashrc whenever you have changed ~/.bashrc).  Also appen your
PYTHONPATH in your ~/.bashrc (or ~/.profile)::

    alias python="~/.local/bin/python"
    PYTHONPATH="~/.local/lib:$PYTHONPATH"
    
    
2. numpy
^^^^^^^^

Directions sourced from `here <https://software.intel.com/en-us/articles/numpy-scipy-with-mkl>`_
Untar and go into directory.::

    cd ~/packages
    tar xzvf numpy-1.8.1.tar.gz
    cd numpy-1.8.1

Add these lines to the site.cfg file to tell numpy where MKL is::

    [mkl]
    library_dirs = /opt/apps/intel/2013.1.039/mkl/lib/intel64
    include_dirs = /opt/apps/intel/2013.1.039/mkl/include
    mkl_libs = mkl_rt
    lapack_libs =

Change the compiler options in the __init__ function of the 
class ``IntelEM64TCCompiler`` in the file: numpy/distutils/intelcompiler.py
Change self.cc_exe to:::

    self.cc_exe = 'icc -O3 -g -fPIC -fp-model strict -fomit-frame-pointer -openmp -xhost'

I skipped the step about changing the fortran compiler because it is not
straightforward and may already be configured adequately.  Now compile and
install. ::


    python setup.py config --compiler=intelem build_clib --compiler=intelem build_ext --compiler=intelem install --user --record files.txt

Explaination of flags: ``--user`` Tells python to install in ~/.local.
``--record files.txt`` Leaves a record of all installed files in case you need
to delte them later ``--compiler=intelem`` Tells python to compile for Intel
64-bit architecture.  Hopefully that all went well. You can optionally add "&>
install.log" to the end of the above command to have the output put in a file
for inspection. You should be able to import it now from any working directory.


3. scipy
^^^^^^^^

Continue directions from numpy source above.


Untar and install.::

    cd ~/packages
    tar xzvf scipy-0.14.0.tar.gz
    cd scipy-0.14.0
    python setup.py config --compiler=intelem --fcompiler=intelem build_clib --compiler=intelem --fcompiler=intelem build_ext --compiler=intelem --fcompiler=intelem install --user --record files.txt


Always good to return to your home directory and try to import a newly
installed package to test if its okay. Also starting the python interpreter
with "-v" will show which directory packages are being loaded from. 
    
4. matplotlib
^^^^^^^^^^^^^

Untar and install.::

    cd ~/packages
    tar xzvf matplotlib-1.3.1.tar.gz
    cd matplotlib-1.3.1
    python setup.py install --user --record files.txt

    
5. ipython 
^^^^^^^^^^

Untar and install.::

    cd ~/packages
    tar xzvf ipython-2.1.0.tar.gz
    cd ipython-2.1.0
    python setup.py install --user --record files.txt


6. pandas
^^^^^^^^^

This is a dependency of MDTraj. Untar and install.::

    cd ~/packages
    tar xzvf pandas-0.13.1.tar.gz
    cd pandas-0.13.1
    python setup.py install --user --record files.txt


7. Cython
^^^^^^^^^


Cython is used to get compiled C-extensions that can drastically improve
bottlenecks. Untar and install. ::

    cd ~/packages
    tar xzvf Cython-0.20.1.tar.gz
    cd Cython-0.20.1
    python setup.py install --user --record files.txt

8. nose
^^^^^^^

Used for unittests. Untar and install.::

    cd ~/packages
    tar xzvf nose-1.3.3.tar.gz
    cd nose-1.3.3
    python setup.py install --user --record files.txt


9. mdtraj 
^^^^^^^^^

Used for manipulating trajectories (e.g. can load ``xtc`` format). Unzip and
install.::

    cd ~/packages
    unzip mdtraj-master.zip
    cd mdtraj-master 
    python setup.py install --user --record files.txt


10. mpi4py
^^^^^^^^^^

Used for running parallel code. When running any program with mpi4py you need to 
load openmpi libraries with:::

    module load openmpi/1.4.4-intel

Untar and install::

    cd ~/packages
    tar xzvf mpi4py-1.3.1.tar.gz
    cd mpi4py-1.3.1
    python setup.py install --user --record files.txt


11. modeller
^^^^^^^^^^^^

Used to make mutant pdbs structures. Untar and install:::

    cd ~/packages


FILL IN DETAILS
