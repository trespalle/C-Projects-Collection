The aim of this code is to estimate the values of the couplings of a Viana-Bray spin glass from equilibrium configurations of it.
A VB spin glass is a spin glass in which the spins are connected through an Erdös-Renyi graph. In our particular case, each coupling has a value of +1
with probability p, and of -1 with probability 1-p.

Given a temperature T, the code inferes the values of the couplings from equilibrium configurations of the system at that temperature, and then computes 
gamma (which is a measurement of the deviation between the estimates and the actual values of the couplings). The value of gamma is averaged over different
realizations, and the whole process is repeated for a set of different values of T. At the end, a gamma vs. T curve is obtained.

To estimate the values of the couplings from the equilibrium configurations, the code finds the values of the couplings that maximize the "pseudo-likelihood 
function" (PLH), which is a function of the sampled configurations (and, of course, of the couplings' estimates). This maximization is done through simulated annealing.

To run the code, the following inputs are needed:

- A file called "input_PLH.txt", whith the values of the relevant parameters:
N (number of spins)
C (number of sampled equilibrium configurations per temperature and per seed)
n_temps (number of temperature values for which gamma is going to be computed)
temp_1, ..., temp_ntemps (values of those temperatures, each of them expressed as a three-digit integer)
n_seeds (number of seeds or independent realizations per temperature)
teff_ini (initial value of the effective temperature used in the simulated annealing)
tau (number of steps of the simulated annealing)
firs_seed (in the simulated annealing calculations, only the configurations corresponding to seeds between first_seed and (first_seed + n_seeds -1) will be used)


- A set of folders named "tXXX", each of which contains equilibrium configurations of the VB spin glass sampled at temperature XXX*0.01
(for example "t125" contains equilibrium configurations sampled at T = 1.25). Each "tXXX" folder contains several binary files named "D_YYY.dat", where
"YYY" is a three-digit number between 000 and n_seeds-1 (that is, if n_seeds = 100, between 000 and 099, for example).
Each "D_YYY.dat" corresponds to an independent realization of the MC simulation from which the equilibrium configurations are sampled (different seed of the RNG), 
and stores C different equilibrium configurations of the spin glass, one per each row (each row is a sequence of N binary values,
+1 or -1, corresponding to "up" and "down" spins, that represents the state of the system). Each "D_YYY.dat" corresponds, as well, to a different realizaition
of the spin glass graph.

- A folder called "GEOMETRY" that contains n_seeds binary files named "GEO_YYY.dat", each of which contains information about the YYY-th realization of the spin 
glass graph. The content of each "GEO_YYY.dat" binary file is: M (the number of bonds), zmax (the maximum connectivity of the graph), neigbhour list row by row.


The output of the code is a file containing the gamma vs. T curve. The file contains one row per each temperature, each row of the form:
temperature, sum_gamma, sum_gamma_2, n_seeds

sum_gamma/n_seeds is the avarage value of gamma for a given temperature. 
((sum_gamma_2/n_seeds) - (sum_gamma/n_seeds)^2)/(n_seeds-1) is our estimate of the error of the average value of gamma at a given temperature.

 

