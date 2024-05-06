***********************  WARNING ***********************

This code does not initialize the Dyson matrix very efficiently. This limits the maximum dimension of the model space to D < 20 (at least on my 16GB laptop). I never needed to go to model spaces larger than D = 10, however, if you do need large model spaces for your research project (or any other reason) please email me at stefano.brolli@unimi.it and we can work together on a more efficient extension of the code that can handle larger model spaces.

********************************************************







****************** TERMS AND CONDITIONS ******************

- The software is provided as it is, without any guarantee or commitment for support. 
 
- You are welcome to use and redistribute this software, as long as it is distributed as a whole and this file is always included and unmodified from its original version. If the software is redistributed with any modification from its original form, this should be clearly stated.

- When using the software for research purposes, I kindly ask you to acknowledge my contributions with the creation of this code.

*********************************************************






******************* WHAT THE CODE DOES **************************

The Richardson Hamiltonian H^{(D)} for D energy levels is defined as 

H^{(D)} = H_0 + H_int, with

H_0 = \xi \sum_{p = 1}^{D} \sum_{s = \uparrow, \downarrow} \left( p - 1 \right) c_{p s}^{\dagger} c_{p s} 

H_int = - \frac{g}{2} \sum_{p, q = 1}^{D} c_{p \uparrow}^{\dagger} c_{p \downarrow}^{\dagger} c_{q \downarrow} c_{q \uparrow}

Please notice the factor 1/2 in front of the interaction, as in some conventions it is absent. We will further assume xi = 1. 
If you need to change it you can edit the file "Main.cpp" and change the variable xi as required.


This code calculates the single-particle Dyson propagator (for closed-shell systems, i.e. an even number of particles here) of the Richardson model in the following approximation schemes: 

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

The source code is found in the folder "source". To install, modify the Makefile to use your C++ compiler and run "make" in the source path. 
The compilation also requires the Eigen 3.4 library. To link the Eigen libraries correctly please modify the "Basis.h" source so that the Eigen Dense and Sparse libraries are properly included.

*********************************************************





******************* CODE USAGE **************************

Program usage: ./RichModSCGF D A g approx start_prop

The approximation scheme (approx) can be ADC2, TDA, ADC3.

The HF propagator is always computed.

The starting propagator (start_prop) can be unperturbed or HF. To use the unperturbed propagator write unpert, anything else will do HF.

To require OpRS, add anything at the end of the line.

Examples:

./RichModSCGF 10 10 0.3 ADC2 unpert         --------> Performs ADC2 + sc0 with an unperturbed reference for D = A = 10, g = 0.3
./RichModSCGF 6 4 0.5 ADC3 HF OpRS          --------> Performs ADC3 + OpRS + sc0  starting with a HF reference for D = 6, A = 4, g = 0.5
./RichModSCGF 8 2 -0.2 TDA HF cats          --------> Performs TDA + OpRS + sc0  starting with a HF reference for D = 8, A = 2, g = -0.2

*********************************************************





******************* CODE OUTPUT **************************

The code outputs on screen the final propagator (calculated with the chosen approximation scheme), the G.S. energy, and the correlation energy. 
It also outputs the same quantities for the HF propagator, regardless of the approximation scheme chosen.

The propagators are displayed as s.p. removal (addition) energy and corresponding spectroscopic amplitude. 

The final propagator is also printed on file in the SpProp/ folder. The hole part is stored in "spprop_backward_sc0.dat" and the particle part is stored in "spprop_forward_sc0.dat". The first element of both files is the number of s.p. energies.

The intermediate propagators calculated during the last sc0 are stored in "spprop_forward_sc0_itr*.dat" and "spprop_backward_sc0_itr*.dat". 
The intermediate propagators generated during the OpRS procedure are not stored.

The HF propagator is stored in SpProp/HF/ as "spprop_backward_HF.dat" and "spprop_forward_HF.dat". 

The final OpRS propagator is stored into SpProp/OpRS/ as "spprop_backward_OpRS.dat" and "spprop_forward_OpRS.dat", with the same notation as the final propagator to distinguish between particles and holes.


The code also outputs the spectral function obtained with the final propagator into SpecFunc/.
A file called "Model_Parameters.dat" is also produced. This file contains in order D, A, g, xi, approximation scheme, starting propagator, and self-consistent scheme. This is useful to know the parameters of the last run and what is contained in the SpProp/ and SpecFunc/ folders.

**********************************************************
