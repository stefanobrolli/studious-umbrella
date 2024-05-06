#include <cmath>
#include <iostream>
#include <iomanip>

#include "Basis.h"

using namespace std;

SpBasis::SpBasis(const unsigned int D, const unsigned int A, double xi): e_0(D) {

	Num_Levels = D;
	Num_Holes = A;

	cout << "The basis is:" << endl;
	cout << endl;
	
	cout << "alpha" << setw(12) << "spin" << setw(15) << "e_0" << endl;
	
	for (int ialpha = 0; ialpha < D; ++ialpha) {
			
		e_0(ialpha) = xi*ialpha;
		
		for (int ispin = 0; ispin < 2; ++ispin) {

			cout << ialpha << setw(15) << 2*ispin - 1 << setw(15) << e_0(ialpha) << endl;
		}
	}
		
	cout << "Number of holes: A = " << Num_Holes << endl;
	cout << "Number of different unperturbed energies: D = " << Num_Levels << endl;
}
	
SpBasis::~SpBasis() {

}


