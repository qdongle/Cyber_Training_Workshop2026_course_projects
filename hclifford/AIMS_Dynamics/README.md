This folder contains all of the files relevant to performing the ab initio multiple spawning dynamics in PySPAWN

There are individual directories for each initial condition

start.py is the input file of the trajectory. It has all of the dynamics set up details such as timestep, legnth of trajectory, and other details of the calculation

start.log is the output file

sim.hdf5, sim.json, sim2.hdf5, and sim2.json are files generated in PySPAWN that store info about the AIMS wave function, energies, forces, position and momenta of each trajectory basis function, and other quantities

geometry.xyz is the xyz coordinates of the initial condition

Slurm.log is the slurm file used to submit the trajectory

hessian.hdf5 is the Hessian used for Wigner sampling of positions and momenta
