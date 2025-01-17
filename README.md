# COMP_SQUARE

This is a compressible Navier-Stokes Solver written in Fortran. It was written under the guidance of Dr. Nagabhushan Rao Vadlamani in his Advanced CFD course taught at IIT Madras. It is based on another high-order solver Dr NRV wrote and distributed as part of the course. Current boundary conditions and initial values are set to solve for a standard Taylor-Green Vortex case for benchmarking.

An explicit 4th-order Runga-Kutta method is used for time marching. The Fluxes are calculated at each of the 4 RK sub-steps and discretized using Explicit (2nd and 4th order) or Compact (4th and 6th order) schemes. Filtering schemes of 4th and 6th order are used to eliminate numerical instabilities. 


