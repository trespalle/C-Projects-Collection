This code runs a continuous-time Monte Carlo simulation of the 2D Ising model in a square lattice with periodic boundary conditions.
The time the system spends in each configuration is a random variable that is sampled without rejection.
(The continuous-time Monte Carlo method is also known as "kinetic Monte Carlo" or "Gillespie" algorithm).

The input variables are specified between lines 47 and 51:
T (temperature)
L (linear size of the 2D lattice)
n_mcs (number of regular, finite MC steps we want our CTMC simulation to last at least)

The output of the program is a file called "evolution_CTMC.txt" which contains three columns:
t (expressed in MC steps), energy of the system, total magnetization of the system

A snapshot of the system after the simulation is generated aswell: a file called "snapshot.txt", where each row is a 
row of the lattice, each "+1" represents a spin pointing "up" and each -1 represents a spin pointing "down".
