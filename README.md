## Simulation: reliable lower bound for QKD protocols key rates
This repo contains Python and Fortran 90 code to perform a reliable lower bound for QKD key rate.
Theory is briefly written in the main files or there's a full documentation in the 'bibliography' directory.
The most useful are
 - https://doi.org/10.22331/q-2018-07-26-77
 - https://doi.org/10.1103/PhysRevResearch.3.013274

# Python scripts
The repo contains a 'brute force' (in the sense that it contains all the functions and operators without calling an external class) script called 'main.py'.
This may be useful to follow all the logic of the simulation.

On the other hand, if one searches a more flexible software will find a Python class called QKD in the file 'src/qkd.py'.
Some examples for the usage of this class can be found in the main directory (simple_BB84.py, main_bb84_4.py, etc.).

For what concerned the SDP solver, I recommend using 'MOSEK' solver (it requires a license).
It is ok to use 'CVXOPT' instead, but it crashes quite often.

# Fortran 90
The simulation is also made in Fortran 90. 
The files are main.f90 and simple_BB84.f90 which requires the modules (in 'src' directory):
 - debugging.f90, where are defined checkpoints and print functions;
 - matrices.f90, useful for matrices calculations (over all the tedious logarithm of a matrix)
 - QKD.f90, containing functions and subroutines calculating entropies and SPDA.
I have used SDPA solver, which is an executable file SPDA.exe which requires an input *.dat file and writes the solution onto a *.out file (download https://sdpa.sourceforge.net/download.html).
SDPA.exe must be in the same folder of the executables.

Still Fortran 90 program has problems with the SDP solution, I'll try to fix them.