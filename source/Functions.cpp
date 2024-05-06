#include "Functions.h"


// Pretty view of complex numbers 

ostream& operator << (ostream& os, const complex <double> z) {
	
	if( abs(imag(z)) < 10E-10 ) {
		
		os << real(z);
	}
	else if( imag(z) < 0 ) {
		
		if ( abs(real(z)) < 10E-10 ) {

			os << "-i" << abs(imag(z));
		}
		else {

			os << real(z) << "-i" << abs(imag(z));
		}
	}
	else {

		if ( abs(real(z)) < 10E-10 ) {

			os << "i" << abs(imag(z));
		}
		else {

			os << real(z) << "+i" << abs(imag(z));
		}
	}
	return os;
}
