The Duplex HCO model is a compartmental model (like SIRS, for example) that studies the propagation of norm violation in a 
population. A detailed explanation and mathematical description of this model can be found in "Duplex_HCO.pdf", which is my Bachelor's 
thesis. In this model, population is divided into 5 different compartments.

By making a few assumptions, one can build a system of 4 nonlinear differential equations for this model. These 4 differential 
equations (that constitute a nonlinear dynamical system) dictate the time evolution of the number of individuals in each compartment.
The 4D dynamical system depends on 6 parameters: alpha, beta, r, p, k1 and k2 (the meaning of these parameters can be found in my Bachelor's thesis).

This code finds all the fixed points of the dynamical system that lie within a certain region of the phase space for several values 
of alpha and fixed values of the rest of the parameters. To do so, the Newton-Raphson algorithm is used.


The parameters beta, r, p, k1 and k2 are definded as macros, as well as n_explorers (see below).

Then, the core of the code is as follows:


1. Alpha is set to some initial value.

2. The positions of the n_explorers "explorers" are initialized in such a way that they are evenly distributed throughout the
region of interest. Each "explorer" can be thought as an "imaginary" point-like particle that will be carried to the nearest fixed
point by the Newton-Raphson algorithm.

3. The Newton-Raphson algorithm is applied to each one of the n_explorers "explorers". Each explorer ends up at a  different fixed point
(some of the explorers end up at fixed points that are out of the region of interest). If n_explorers is sufficiently large, every relevant
fixed point of the system inside the region of interest is reached by at least one explorer.

4. The final positions of the explorers are compared, and the ones that are different (that correspond to different fixed points)
are stored in the matrix "descubrimientos". Each row of the "descubrimientos" matrix is the 4D position of one of the discovered 
fixed points.

5. The stability of each discovered fixed point is computed, and stored in the fifth column of "descubrimientos" matrix.
(stability = -1 if the fixed point is stable, and +1 otherwise).

5. The discovered fixed points are appended to a file called "Raphson_lines.txt" in the format: alpha, x, y, z, w, stability
((x, y, z, w) is the 4D position of the fixed point).

6. The value of alpha is increased.

7. Back to step 2.

The program ends when alpha = 1.0.


The output of the program is the file "Raphson_lines.txt", that contains all the fixed points inside the region of interest
for all the explored values of alpha between 0 and 1. Each row of "Raphson_lines.txt" is of the form (alpha, x, y, z, w, stability),
where (x, y, z, w) is the 4D position of the fixed point and stability = -1 if the fixed point is stable (and +1 otherwise).


 