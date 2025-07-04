# LES3D_TurbulentDiffusionFlame_MATLAB
This repository provides a MATLAB code for simulating turbulent diffusion flames at low Mach numbers using Large Eddy Simulation (LES). The solver employs the finite difference method for spatial discretization and a projection-type algorithm for time integration. A companion article is under submission to a peer-reviewed journal.

# Running the Simulation
To run the simulation in MATLAB, ensure that all the following auxiliary functions are available in the working directory. Each one is responsible for a specific physical or numerical computation:

# Core Numerical Routines
convective2.m
Computes convective terms using second-order central finite differences.

convective4.m
Computes convective terms using fourth-order central finite differences.

pressure2.m
Calculates the pressure gradient using second-order finite differences.

pressure4.m
Calculates the pressure gradient using fourth-order finite differences.

diffusive2.m
Computes viscous/diffusive terms using second-order finite differences.

diffusive4.m
Computes viscous/diffusive terms using fourth-order finite differences.

poisson2.m
Solves the Poisson equation for pressure correction using a second-order discretization.

poisson4.m
Solves the Poisson equation for pressure correction using a fourth-order discretization.

tempdens.m
Determines temperature and density fields based on the mixture fraction, using thermodynamic relations and flame assumptions.

animation.m
Generates an animated .gif showing the temporal evolution of key flow variables (velocity, temperature, mixture fraction), aiding visualization and analysis.
