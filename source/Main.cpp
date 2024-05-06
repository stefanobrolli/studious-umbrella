// We focus on the spin up propagator, by simmetry the spin down one is the same.

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>


#include "Basis.h"
#include "Props.h"
#include "Functions.h"
#include "SelfEn.h"


using namespace std;

namespace fs = std::filesystem;

int main(int argc, char **argv) {
	
	unsigned int D, A;  
	double g; 
	double xi = 1.; // We set the distance between unperturbed levels to be 1. CHANGE IT IF NEEDED 
	string approx;
	string start_prop;

	// Creating the necessary directories

	mkdir("./SpecFunc", 0777);
	mkdir("./SpProp", 0777);
	mkdir("./SpProp/HF", 0777);
	mkdir("./SpProp/OpRS", 0777);

	if( argc == 6 || argc == 7 ) {
	
		D = stoi(argv[1]);
		A = stoi(argv[2]);
		g = stod(argv[3]);
		approx = argv[4];
		start_prop = argv[5];
	}
	else {

		cerr << "Program usage: ./RichModSCGF D A g approx start_prop" << endl;
		cerr << "The approximation scheme (approx) can be ADC2, TDA, ADC3" << endl;
		cerr << "The starting propagator (start_prop) can be unperturbed or HF. To use the unperturbed propagator write unpert, anything else will do HF." << endl;
		cerr << "To require OpRS add anything at the end of the line." << endl;
		cerr << "Examples:" << endl;
		cerr << "./RichModSCGF 10 10 0.3 ADC2 unpert" << endl;
		cerr << "./RichModSCGF 10 10 0.3 ADC3 HF OpRS" << endl;
		exit(4);
	}

	cout << setprecision(10);

	cout << "We set xi = " << xi << " and g = " << g << endl;
	
	// Printing this info on a file
	
	string paramfile = "Model_Parameters.dat";		
    
	ofstream fout(paramfile);		
				
	fout << D << endl;
	fout << A << endl;
	fout << g << endl;
	fout << xi << endl;
	fout << approx << endl;
	fout << start_prop << endl;

	if (argc == 6) {

		fout << "sc0" << endl;
	}
	else if (argc == 7) {

		fout << "OpRS + sc0" << endl;
	}
 	
 	fout.clear();
 	fout.close();
	
	SpBasis SpBas(D, A, xi);

	SpProp G0;
	G0.Build_G0(SpBas);

	
	SelfEn SE_HF(SpBas, G0, g);

	SE_HF.BuildStatic(G0);
	SE_HF.Tadpole(SpBas, "./SpProp/HF/", "HF");

	
	SpProp G;
	G.ChargeNewProp(SpBas, "./SpProp/HF/", "HF");

	double EHF = G.MGK(SpBas);

	cout << "The Hartree-Fock energy is: " << EHF << endl;
	cout << endl;

	cout << "The Hartree-Fock propagator is:" << endl;

	G.Print();

	double EGS = 0.;
	double EGS_old = 1.;

	if (argc == 6) {

		if (start_prop == "unpert") {

			G.Build_G0(SpBas);
		}

		//G.Print();

		SelfEn SE(SpBas, G, g);
		SE.BuildStatic(G);
	
		if (approx == "ADC2") {

			SE.BuildADC2(G);	
		}
		else if (approx == "TDA") {

			SE.BuildTDA(G);	
		}
		else if (approx == "ADC3") {

			SE.BuildADC3(G);	
		}
		else {

			cerr << "Approximation not known..." << endl;
			exit(4);
		}

		SE.sc0(SpBas, G);
		G.GenerateOpRS();
		G.ChargeNewProp(SpBas, "./SpProp/OpRS/", "OpRS");
	}	
	else if (argc == 7) {

		int count = 0;

		if (start_prop == "unpert") {

			G.Build_G0(SpBas);
		}

		while (abs(EGS_old - EGS) > 1E-7) {

			cout << "OpRS ITERATION #" << ++count << ":" << endl;
			
			SelfEn SE(SpBas, G, g);
			SE.BuildStatic(G);

			if (approx == "ADC2") {

				SE.BuildADC2(G);	
			}
			else if (approx == "TDA") {

				SE.BuildTDA(G);	
			}
			else if (approx == "ADC3") {

				SE.BuildADC3(G);	
			}
			else {

				cerr << "Approximation not known..." << endl;
				exit(4);
			}
	
			SE.sc0(SpBas, G);

			G.ChargeNewProp(SpBas, "./SpProp/", "sc0");
			//G.Print();
			EGS_old = EGS;
			EGS = G.MGK(SpBas);
		
			G.GenerateOpRS();
			G.ChargeNewProp(SpBas, "./SpProp/OpRS/", "OpRS");
			//G.Print();

			cout << "E_GS = " << EGS << ",  Correlation energy = " << EGS - EHF << endl;
			cout << "Correlation energy captured by OpRS = " << G.MGK(SpBas) - EHF  << endl;

			cout << endl;
			cout << endl;
		}
	}

	double E_OpRS = G.MGK(SpBas);

	cout << "The final propagator is:" << endl;

	G.ChargeNewProp(SpBas, "./SpProp/", "sc0");
	G.Print();
	EGS = G.MGK(SpBas);
	cout << "The final GS energy is E_GS = " << EGS << ".  Correlation energy = " << EGS - EHF << endl;
	cout << "Correlation energy captured by OpRS = " << E_OpRS - EHF  << endl;

	G.SpecFunc();
	
	return 0;
}
