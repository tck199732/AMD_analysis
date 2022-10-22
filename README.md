# amd_master


**Background**
This repository contains fortran codes for simulating heavy ion collision in the framework of 'Antisymmetrized Molecular Dynamics (AMD)' (references) and C++ code for analyzing the output of the simulations. Events and particles are selected and compared to experimental data from NSCL E15190.

**Table of contents**
1. [Installation](#1-installation)
1. [Structure of the repository](#2-structure-of-the-repository)
1. [Testing framework](#3-testing-framework)
1. [PEP 8 style guide](#4-pep-8-style-guide)

## 1. Installation


## 2. Structure of the repository
This repository is developed in Python 3.8 and C++17 (`-std=c++17` in GCC 9.4.0).
- [**`simulation/`**](simulation/): Fortran code for HIC simulations. 
- [**`database/`**](database/): For storing all the experimental data and parameters.
- [**`src/`**](src/): C++ source files for analyzing output of AMD simulation.
- [**`analysis/`**](src/): C++ source files for analyzing output of AMD simulation.
