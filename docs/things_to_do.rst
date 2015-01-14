Things To Do and Keep in Mind
=============================

Things To Keep in Mind
----------------------

Our overall coding goals (in roughly decreasing importance):

1. Write code that is easy to understand and extend: style, design, and documentation.

  - Please force yourself to adhere to `PEP 8 style guide. <http://legacy.python.org/dev/peps/pep-0008>`_
  - Standardizing our coding style promotes readability. 
  - 20 useful adages: `Zen of Python <http://legacy.python.org/dev/peps/pep-0020/>`_
  - Program design and development is nontrivial, ideally we would follow the 
    new school philosophy of `"agile development" <http://en.wikipedia.org/wiki/Agile_software_development>`_
    where development process involves many small iterations. 
  - We should aim to add more features over time while simultaneously reducing our number of source code lines.

2. Write code that works: unit testing.
  - Code is broken until tested.

3. Eventually write code that is efficient with time and memory, but only after #1 & #2.

Things To Do
------------

Things to be improved, tested, or implemented for the first time.

General Objectives
^^^^^^^^^^^^^^^^^^

1. Generally, we want to minimize our dependence on outside packages.

2. Descriptive and concise docstrings in all modules. Decide on essentials
   for docstrings.

3. Figure out how to grab source code docstrings into Sphinx documentation.

4. Write a simple example document.

Models
^^^^^^

1. Submodules for bonded_potentials & pairwise_potentials to store
   contact interaction functions:
    1. Bonded potential utility library. DONE
    2. Pairwise potentials library. With dimensionless 
       pairwise potential types:    DONE
        - LJ1210
        - LJ1210rep
        - Bare Gaussian
        - Gaussian + hard wall
        - Double Gaussian + hard wall
        - Others? (Desolvation barriers, etc.)
    3. Consider allowing for interactions made up of a small set of 
       basis functions.
    4. Determine a quick way to evaluate the energies of contacts. 
    
2. Define the potential energy as the sum of "interaction strengths"
   multiplied by corresponding "interaction potentials". Then the list
   of non-redundant interaction strengths are the "model parameters".
   Want to generalize parameter fitting to act on model parameters.  DONE
   
3. Code disulfides as a bonded interaction with harmonic bond, angle, etc. DONE
   
4. Decide on design that allows for multiple chains in the topology.

5. Deprecated distinquishing "Tf" from "Mut" instead just use iteration. DONE

6. Have bonded interactions set in separate programs that deal with bead representations.
