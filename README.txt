***********************  WARNING ***********************

This code does not initialise the Dyson matrix very efficiently. This limits the maximum dimension of the model space to D < 20 (at least on my 16GB
laptop). I never had the need to go to model spaces larger than D = 10, however if you do need large model spaces for your research project 
(or any other reason) please email me at stefano.brolli@unimi.it and we can work together on a more efficient extension of the code that can handle
larger model spaces.

********************************************************







****************** TERMS AND CONDITION ******************

- The software is provided as it is, without any guarantee or commitment for support. 
 
- You are welcome to use and redistribute this software, as long as it is distributed as a whole and this file is always included and unmodified 
from its original version. If the software is redistributed with any modification from its original form, this should be clearly stated.

- When using the software for research purposes, I kindly ask you to acknowledge my contributions with the creation of this code.

*********************************************************






******************* WHAT THE CODE DOES **************************

The Richardson Hamiltonian H^{(D)} for D energy levels is defined as 

H^{(D)} = H_0 + H_int, with

H_0 = \xi \sum_{p = 1}^{D} \sum_{s = \uparrow, \downarrow} \left( p - 1 \right) c_{p s}^{\dagger} c_{p s} 

H_int = - \frac{g}{2} \sum_{p, q = 1}^{D} c_{p \uparrow}^{\dagger} c_{p \downarrow}^{\dagger} c_{q \downarrow} c_{q \uparrow}

Please notice the factor 1/2 in front of the interaction, as in some conventions it is absent. We will further assume xi = 1. 
If you need to change it you can edit the file "Main.cpp" and change the variable xi as required.


This code calculates the single-particle Dyson propagator (for closed shell systems, i.e. an even number of particles here) of the Richardson model 
in the following approximation schemes: 

1) Hartree-Fock (HF)
2) ADC(2) (second order calculation of the irreducible self-energy)
3) 2p-1h TDA
4) ADC(3)

The self-consistent scheme can be chosen between sc0 (with either unperturbed or HF reference propagator) and OpRS (optimized reference state) + sc0.
For a detailed explanation of sc0 and OpRS please see: 

sc0: https://link.aps.org/doi/10.1103/PhysRevC.89.024323
OpRS: https://journals.aps.org/prc/abstract/10.1103/PhysRevC.105.044330


*****************************************************************







****************** INSTALLATION *************************

The source code is found in the folder "source". To install modify the Makefile to use your C++ compiler and run "make" in source path. 
The compilation also requires the Eigen 3.4 library. To link the Eigen libraries correctly please modify the "Basis.h" source so that the Eigen Dense 
and Sparse libraries are properly included. 

*********************************************************





******************* CODE USAGE **************************

Program usage: ./RichModSCGF D A g approx start_prop

The approximation scheme (approx) can be ADC2, TDA, ADC3

The starting propagator (start_prop) can be unperturbed or HF. To use the unperturbed propagator write unpert, anything else will do HF.

To require OpRS add anything at the end of the line.

Examples:

./RichModSCGF 10 10 0.3 ADC2 unpert         --------> Performs ADC2 + sc0 with an unperturbed reference for D = A = 10, g = 0.3
./RichModSCGF 6 4 0.5 ADC3 HF OpRS          --------> Performs ADC3 + OpRS + sc0  starting with a HF reference for D = 6, A = 4, g = 0.5
./RichModSCGF 8 2 -0.2 TDA HF cats          --------> Performs TDA + OpRS + sc0  starting with a HF reference for D = 8, A = 2, g = -0.2

*********************************************************