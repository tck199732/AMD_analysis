# amd_analysis


**Background**
This repository contains codes for analyzing the output of the 'Antisymmetrized Molecular Dynamics (AMD)' simulations. Events and particles are selected and compared to experimental data from NSCL E15190.

**Table of contents**
1. [Installation](#1-installation)
2. [Structure of the repository](#2-structure-of-the-repository)
3. [Procedures](#2-procedures-of-the-analysis)

## 1. Installation
```bash
git clone https://github.com/tck199732/amd_analysis.git
```

## 2. Structure of the repository
This repository is developed in Python 3.8.10 and C++17 (`-std=c++17` in GCC 9.4.0).
- [**`database/`**](database/): For storing all the experimental data and parameters.
- [**`src/`**](src/): C++ source files for analyzing output of AMD simulation.
- [**`bin/`**](src/): C++ main program for analyzing output of AMD simulation.
- [**`pyamd/`**](src/): python source files for analyzing output from C++ main program
- [**`build.py`**](build.py): For setting up conda environment with ROOT ver 6.26.06


## 3. Procedures
The starting point of the analysis is the raw output file from AMD simulation, i.e. table21.dat and table3.dat. In case we want to analyze the spacetime of last interation, we would also need the hist_coll.dat and amdgid.dat which contains info of collision history and the correspondance between nucleon global ID and primary particle. 

- Convert table21.dat and table3.dat to table21.root and table3.root by using the main program amd2root.cpp.
```bash
cd ${project_dir}/bin
g++ amd2root.cpp -o amd2root -I`root-config --libs --glibs --cflags` -I../src
./amd2root {reaction} {mode}, {path_input}, {path_output}
```
where reaction refers to the reaction system such as 'Ca40Ni58E140', mode refers to analysis mode ('21', '3', '21t'), etc

- It is easy to write a script for analysis for pure simulation without experimental constraint. To compare AMD result with experiment, one needs to filter the events using ExpFilter program. For e15190, run 
```bash
cd {project_dir}/bin
g++ ExpFilter.cpp -o ExpFilter -I`root-config --libs --glibs --cflags` -I../src
./ExpFilter {reaction} {mode}, {data_dir}, {path_list}, {path_output}
```
where data_dir refers to directory of the data in root format, path_list refers to the path of a list which contains all files included in the analysis.

- to be continued.

