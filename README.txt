"polymer_evolution.c" simulates the evolution of a polymer in 3D in contact with a thermal bath (canonical ensemble) by integration
of the Langevin equation. The Langevin equation is integrated using the Runge-Kutta 4 method.
The polymer is modeled as a chain of N_part point-like monomers connected by springs (harmonic oscillators of 
constant k_e).

The values of the relevant parameters are specified in lines 8 to 15:
N_parts (number of monomers), m (mass of each monomer), T (temperature of the thermal bath), b (equilibrium length of each spring),
gamma (each harmonic oscillator has an elongation of gamma*b and zero velocity at the beginning of the simulation, gamma is a 
fraction between 0 and 1), k_e (spring constant).

Some other relevant parameters are initialized in lines 69 to 71:
eta (the damping parameter that appears in the Langevin equation), h (the time step), tiempo_fisico (the duration of the simulation
in "physical" units of time).

The output of "polymer_evolution.c" is:

- A file called "polymer_evolution.txt" with the time evolution of the main observables of the system. This file has one row per 
time step, and 5 columns: time, distance between monomers 0 and 1, potential energy per monomer, kinetic energy per monomer, 
radius of gyration-

- A file called "final_picture.txt" with the positions of each monomer at the end of the simulation (one row per monomer, each row 
of the form x, y, z).





