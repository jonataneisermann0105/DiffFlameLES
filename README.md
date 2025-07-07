# DiffFlameLES

This repository provides a MATLAB code for simulating turbulent diffusion flames at low Mach numbers using Large Eddy Simulation (LES). The solver employs the finite difference method for spatial discretization and a projection-type algorithm for temporal integration. Subgrid-scale turbulence is modeled using the classical Smagorinsky approach. A companion article is currently under submission to a peer-reviewed journal.

## Running the Simulation

To run the simulation in MATLAB, execute the [`simulation_LES3D.m`](./simulation_LES3D.m) script. Make sure all auxiliary functions listed below are available in the working directory, as each one is responsible for a specific numerical computation or post-processing task:

- [`animation.m`](./animation.m): generates an animated `.gif` showing the temporal evolution of key flow variables (velocity magnitude, pressure, mixture fraction, temperature), aiding visualization and analysis;
- [`convective2.m`](./convective2.m): computes convective terms using second-order central finite differences;  
- [`convective4.m`](./convective4.m): computes convective terms using fourth-order central finite differences;    
- [`diffusive2.m`](./diffusive2.m): computes viscous/diffusive terms using second-order finite differences;  
- [`diffusive4.m`](./diffusive4.m): computes viscous/diffusive terms using fourth-order finite differences;  
- [`poisson2.m`](./poisson2.m): solves the Poisson equation for pressure correction using second-order discretization;  
- [`poisson4.m`](./poisson4.m): solves the Poisson equation for pressure correction using fourth-order discretization;
- [`pressure2.m`](./pressure2.m): calculates the pressure gradient using second-order finite differences;  
- [`pressure4.m`](./pressure4.m): calculates the pressure gradient using fourth-order finite differences;  
- [`tempdens.m`](./tempdens.m): computes temperature and density fields based on the mixture fraction and thermodynamic principles.  

