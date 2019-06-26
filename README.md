# Volume of Fluid Phase Coupling (VOF_PC) UDF Library (vof-phase-coupling-library) - for the coupling between Ansys Fluent and ANSYS Mechanical APDL

## Introduction

This is a small UDF (C) library, which can be used to couple between Ansys Fluent and Ansys Mechanical APDL for Multiphase Volume of Fluid (VoF) calculations. It has been written for the simulation of an alternating current (50 Hz) electroslag-remelting (ESR) process, as part of the DFG Projekt (309143753). In this process current is used as primary energy to melt a metal electrode due to the resistance heating of a slag. As a results metal droplets are falling through the slag layer and collecting bellow, the droplets are arbitrarily changing the current pattern and lorentz forces inside the slag region. To simulate this phenomena one have to resolve both the electromagnetic fluxes and the fluid flow.

Since most low frequency electromagnetic simulation software, for various reasons, is implemented in finite element methods and decent fluid flow solvers are usually implemented in finite volume methods, due to the beneficial conservative properties of the method, one usually have to couple between these methods. 

The application of these library is based on the use of Ansys Fluent and Ansys Mechanical APDL (for EMAG). It is a file based coupling method, which is not optimal, but cannot be handled much better without direct access to MPI functions from APDL and UDF code itself, which is not documented by ANSYS in the required form.
The in this library implemented methods are focused on element/cell based property, respectively source term, changes. Topological based coupling implementations or approaches, which require capturing the interfaces between the phases and the meshing of these disjunct phases, as separate volumes, for the solution within the EMAG software, are by nature more accurate in most FEM implementations, but require more computational effort and are furthermore a lot more difficult to implement in a general/robust way.

**Warning: To make this library to work for your specific case(s), you will probably have to make several modifications to this model, so generally it can not be seen as a "ready to use" library, more as a library template, which can be utilized for your needs.**

*If you are generally interested in coupling methods you may also have a look at [EOF-Library](https://github.com/jvencels/EOF-Library) for OpenFOAM and Elmer, or for more general use in ANSYS you may have a look at the documented system coupling capabilities or the [MPCCI Library](https://www.mpcci.de/)* 

## Library overview
Short overview about limitations, prerequisites and functionality.

### Current limitations
  * only one cell zone area can be coupled in current version of the library
  * only nearest neighbor (NN) cell coupling
  * current implementation assumes fixed meshes during calculation
  * Coupled regions must be geometrical identical (size, (units) and location) for the mapping to make sense, there is no translation or rescaling functionality implemented jet.
  * To exit Fluent during coupling iterations you have to cancel calculation during the iterations of a timestep, or kill the fluent process.

### Prerequisites
  * Multiphase (VOF) Fluent case (parallel 2D or 3D)
  * related Case in form of a Ansys Mechanical APDL script
  * you have to be sure, which quantities shall be exchanged, current examples shows:
    * VOF - Fluent to ANSYS
    * Lorentz force and Joule heating sources from ANSYS to Fluent


### Implemented Coupling Methods
Here you will find a short list of currently implemented coupling methods:

#### Nearest Neighbor Coupling (fixed grids)
This coupling method assumes fixed grids which are not changed between the individual computations. The meshes are mapped to each other via nearest neighbor method (*which is currently implemented as brute force algorithm, since computational efficiency is irrelvant for fixed grids, as the mapping is only done once*).

*Functionality is mainly implemented in vof_pc_nn_coupling.c*

#### Todo ...


## Using the Library
This is a short description of how to use the library as it is by now for the ESR process. So to couple Ansys Mechanical APDL with one cell zone in ANSYS Fluent. For this case the exchanged properties for each timestep are:
  * volume fractions to ANSYS APDL
  * Lorentz forces and Joule heating source terms to Ansys Fluent

**Warning: The following text is no complete instruction. If you want to adopt these library you probably have to understand a lot of the code by yourself anyhow...**

### Modifications to ANSYS APDL Script
Material properties, boundary conditions or the mesh etc. can be changed quite easily.
 If only a specific region of the EMAG simulation should be coupled with the Fluent cell zone you can select these inside at the marked area inside the APDL script region termed as "DIMENSIONING OF NECESSARY VARIABLES/ARRAYS/MATRICES". 

Variable | Description
--- | --- 
`XC_PATH(1)` | has to be adapted to current exchange folder path

**Info: If you want to build your own case a basic template file is given with apdl_example.ans"**

### Modifications ANSYS Fluent UDF Files
At least 3-4 (2D/3D) UDMI's at cell center location will be needed. The main changes can be set within "vof_pc_main.h", once done recompile and load library.

#### vof_pc_main.h
Variable | Description
--- | --- 
`FLUID_ID` | has to be set to cell zone id which shall be coupled with ANSYS APDL 
 `_XC_FOLDER_PATH_` | has to be adapted to current exchange folder path


### Initializing and running the coupling
To initialise the coupling you have to setup your single cases and make the necessary modifications if desired. Bellow you will find the instructions for the selected method of choice *(currently only one available...)*

#### Nearest Neighbor Coupling (fixed grids)
  * start Ansys Fluent and load case and data
  * compile and load this udf library
  * to reset the coupling run the define-on-demand function "ResetSyncState_oD" in Fluent
  * run ANSYS APDL script for example via the powershell file "run-ansys.ps1"
  * if Ansys element coordinates have been written to XC folder and sync.txt status is equal to 1 run define-on-demand function "initCouplingWithAnsysCoords_oD" in Fluent
  * make sure that source terms (Lorentz forces and Joule heating) are enabled inside the coupled cell zone
  * to start the coupling, make sure that "strictCoupling_aE" as execute-at-end function in Fluent is enabled and start the Fluent iterations


## Todo
  * improvement of (code) documentary
  * further code refurbishment and enhancement of code modularity
  * improvement of current limitations
  * make NN faster (i.e. kdtree search)
  * make use of better (more accurate) coupling methods than NN
  * make coupling function for variating grids -> fast mapping necessary



