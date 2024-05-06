#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>

#include "Basis.h"
#include "Props.h"
#include "Functions.h"

using namespace std;


SpProp::SpProp() {

}
	
SpProp::~SpProp() {
	
}
	

SpProp::SpProp(const SpProp& G) {

	Y_kp = G.Y_kp;
	X_np = G.X_np;
	
	esp_k = G.esp_k;
	esp_n = G.esp_n; 
}


SpProp& SpProp::operator = (const SpProp& G) {

	Y_kp = G.Y_kp;
	X_np = G.X_np;
	
	esp_k = G.esp_k;
	esp_n = G.esp_n; 
	
	return *this;
}


void SpProp::Build_G0(const SpBasis& SpBas) {	

	int num_pairs = SpBas.Num_Holes/2;
	
	Y_kp.setZero(num_pairs, SpBas.Num_Levels);
	X_np.setZero(SpBas.Num_Levels - num_pairs, SpBas.Num_Levels);

	esp_k.setZero(num_pairs);
	esp_n.setZero(SpBas.Num_Levels - num_pairs);
	
	// Building the unperturbed propagator
	
	for (int itr_pole = 0; itr_pole < SpBas.Num_Levels; ++itr_pole) {	
			
		if (itr_pole < num_pairs) {

			esp_k(itr_pole) = SpBas.e_0(itr_pole);
			Y_kp(itr_pole, itr_pole) = 1.;
		}
		else {

			esp_n(itr_pole - num_pairs) = SpBas.e_0(itr_pole);
			X_np(itr_pole - num_pairs, itr_pole) = 1.;	
		}
	}
}


void SpProp::Print() {

	// Print backward part of the propagator in partial waves
	
	cout<<endl;
	cout<<endl;
	cout<<"Backward part of the propagator:"<<endl;
	cout<<endl;


	for (unsigned int itr_h = 0; itr_h < Y_kp.rows(); ++itr_h) {

		cout << "e_k = " << esp_k(itr_h) << " :" << endl;

		cout<<"Y_kp = ( ";
		
		for (unsigned int itr_alpha = 0; itr_alpha < Y_kp.cols(); ++itr_alpha) {
		
			cout << Y_kp(itr_h, itr_alpha) << " ";
		}
						
		cout<<")";
		cout<<endl;
		cout<<endl;
	}

	
	// Print forward part of the propagator in partial waves


	cout<<endl;
	cout<<endl;
	cout<<"Forward part of the propagator:"<<endl;
	cout<<endl;


	for (unsigned int itr_p = 0; itr_p < X_np.rows(); ++itr_p) {

		cout << "e_n = " << esp_n(itr_p) << endl;

		cout<<"X_np = ( ";
		
		for (unsigned int itr_alpha = 0; itr_alpha < X_np.cols(); ++itr_alpha) {
		
			cout << X_np(itr_p, itr_alpha) << " ";
		}
						
		cout<<")";
		cout<<endl;
		cout<<endl;
	}
}


void SpProp::ChargeNewProp(const SpBasis& SpBas, const string path, const string id) {

	unsigned int num_poles;
	double esp_i;
	double Z_i;
			
	ifstream fin_back(path + "spprop_backward_" + id + ".dat");

	if(!fin_back) {

		cerr << "Couldn't load the backward part of the new propagator..."<< endl;
		cerr << path + "spprop_backward_" + id + ".dat" << endl;
		exit(7);
	}

	fin_back >> num_poles;

	Y_kp.setZero(num_poles, SpBas.Num_Levels);
	esp_k.setZero(num_poles);

	for(int itr_pol = 0; itr_pol < Y_kp.rows(); ++itr_pol) {

		fin_back >> esp_i;
		esp_k(itr_pol) = esp_i;

		for (int itr_alpha = 0; itr_alpha < Y_kp.cols(); ++itr_alpha) {

			fin_back >> Z_i;
			Y_kp(itr_pol, itr_alpha) = Z_i;
		}
	}
				
	fin_back.clear();
	fin_back.close();

				
	ifstream fin_fwd(path + "spprop_forward_" + id + ".dat");

	if(!fin_fwd) {

		cerr << "Couldn't load the forward part of the new propagator..." << endl;
		cerr << path + "spprop_forward_" + id + ".dat" << endl;
		exit(8);
	}

	fin_fwd >> num_poles;

	X_np.setZero(num_poles, SpBas.Num_Levels);
	esp_n.setZero(num_poles);

	for(int itr_pol = 0; itr_pol < X_np.rows(); ++itr_pol) {

		fin_fwd >> esp_i;
		esp_n(itr_pol) = esp_i;

		for (int itr_alpha = 0; itr_alpha < X_np.cols(); ++itr_alpha) {

			fin_fwd >> Z_i;
			X_np(itr_pol, itr_alpha) = Z_i;
		}
	}

	fin_fwd.clear();
	fin_fwd.close();
	
	//cout<<"New sp propagator built!"<<endl;		
}



void SpProp::PrintFile(const SpBasis& SpBas, const string path, const string id) {
		
	ofstream fout_back(path + "spprop_backward_" + id + ".dat");
	fout_back.precision(15);

	fout_back << Y_kp.rows() << endl;

	for(int itr_pol = 0; itr_pol < Y_kp.rows(); ++itr_pol) {

		fout_back << esp_k(itr_pol) << endl;

		for (int itr_alpha = 0; itr_alpha < Y_kp.cols(); ++itr_alpha) {

			fout_back << Y_kp(itr_pol, itr_alpha) << " ";
		}

		fout_back << endl;
	}
				
	fout_back.clear();
	fout_back.close();


				
	ofstream fout_fwd(path + "/spprop_forward_" + id + ".dat");
	fout_fwd.precision(15);

	fout_fwd << X_np.rows() << endl;

	for(int itr_pol = 0; itr_pol < X_np.rows(); ++itr_pol) {

		fout_fwd << esp_n(itr_pol) << endl;

		for (int itr_alpha = 0; itr_alpha < X_np.cols(); ++itr_alpha) {

			fout_fwd << X_np(itr_pol, itr_alpha) << " ";
		}

		fout_fwd << endl;
	}
				
	fout_fwd.clear();
	fout_fwd.close();
}



void SpProp::SpecFunc() {

	double spec = 0.;
	
	// Print holes

	ofstream fout_h("./SpecFunc/SpectralFunction_holes.dat");
	fout_h.precision(15);

	for(unsigned int itr_h = 0; itr_h < Y_kp.rows(); ++itr_h) {

		for(int itr_alpha = 0; itr_alpha < Y_kp.cols(); ++itr_alpha) {

			spec += pow( Y_kp(itr_h, itr_alpha) , 2);
   		}

		fout_h << esp_k(itr_h) << "  " << spec << endl;
		spec = 0.;
	}

	fout_h.clear();
	fout_h.close();


	// Print particles

	ofstream fout_p("./SpecFunc/SpectralFunction_particles.dat");
	fout_p.precision(15);

	for(unsigned int itr_p = 0; itr_p < X_np.rows(); ++itr_p) {

		for(int itr_alpha = 0; itr_alpha < X_np.cols(); ++itr_alpha) {

			spec += pow(X_np(itr_p, itr_alpha), 2);
   		}

		fout_p << esp_n(itr_p) << "  " << spec << endl;
		spec = 0.;
	}

	fout_p.clear();
	fout_p.close();
}


double SpProp::MGK(const SpBasis& SpBas) const {
	
	double En = 0.;
	
	for(int itr_h = 0; itr_h < Y_kp.rows(); ++itr_h) {
		
		for(int itr_alpha = 0; itr_alpha < SpBas.Num_Levels; ++itr_alpha) {
			
			En += ( esp_k(itr_h) + SpBas.e_0(itr_alpha) )*pow( Y_kp(itr_h, itr_alpha), 2.);
		}
	}

	return En;
}


void SpProp::GenerateOpRS() {

	unsigned int D = X_np.cols();
	double E_F = ( esp_k( esp_k.size() - 1 ) + esp_n(0) )/2.; // esp_k.size() - 1

	DenseVec esp_k_OpRS(D);
	DenseVec esp_n_OpRS(D);

	DenseMat Residues_k_OpRS(D, D);
	DenseMat Residues_n_OpRS(D, D);

	Residues_k_OpRS.setZero();
	Residues_n_OpRS.setZero();


	for (unsigned int alpha = 0; alpha < D; ++alpha) {
		
		// Hole Part
		
		double M_0 = 0.;
		double M_1 = 0.;

		for (int k = 0; k < Y_kp.rows(); ++k) {

			M_0 += pow(Y_kp(k, alpha), 2);
			M_1 += pow(Y_kp(k, alpha), 2)/(E_F - esp_k(k));
		}

		esp_k_OpRS(alpha) = E_F - M_0/M_1;
		Residues_k_OpRS(alpha, alpha) = sqrt(M_0);


		// Particle part

		M_0 = 0.;
		M_1 = 0.;

		for (int n = 0; n < X_np.rows(); ++n) {

			M_0 += pow(X_np(n, alpha), 2);
			M_1 += pow(X_np(n, alpha), 2)/(E_F - esp_n(n));
		}

		esp_n_OpRS(alpha) = E_F - M_0/M_1;
		Residues_n_OpRS(alpha, alpha) = sqrt(M_0);
	}


	// Order the hole poles 

	DenseVec sorted_holes_poles(esp_k_OpRS);
	DenseMat sorted_holes_residues(D, D);

	sort(sorted_holes_poles.data(), sorted_holes_poles.data() + sorted_holes_poles.size());

	for (int itr = 0; itr < sorted_holes_poles.size(); ++itr) {
        
		Eigen::Index index;
		(esp_k_OpRS.array() - sorted_holes_poles(itr)).abs().minCoeff(&index); // Find index of closest element
        sorted_holes_residues.row(itr) = Residues_k_OpRS.row(index); // Rearrange rows
    }
	
	// Order the particle poles 

	DenseVec sorted_part_poles(esp_n_OpRS);
	DenseMat sorted_part_residues(Residues_n_OpRS.rows(), D);

	sort(sorted_part_poles.data(), sorted_part_poles.data() + sorted_part_poles.size());

	for (int itr = 0; itr < sorted_part_poles.size(); ++itr) {
        
		Eigen::Index index;
		(esp_n_OpRS.array() - sorted_part_poles(itr)).abs().minCoeff(&index); // Find index of closest element
        sorted_part_residues.row(itr) = Residues_n_OpRS.row(index); // Rearrange rows
    }

	// Printing it
	
	//cout << "Printing OpRS propagator..." << endl;

	ofstream fout_h("./SpProp/OpRS/spprop_backward_OpRS.dat");
	fout_h.precision(15);

	fout_h << D << endl;

	for (unsigned int k = 0; k < D; ++k) {
			
		fout_h << sorted_holes_poles(k) << endl;

		for (unsigned int itr_comp = 0; itr_comp < D; ++itr_comp) {

			fout_h << sorted_holes_residues(k, itr_comp) << "  ";
		}

		fout_h << endl;
	}

	fout_h.clear();
	fout_h.close();



	ofstream fout_p("./SpProp/OpRS/spprop_forward_OpRS.dat");
	fout_p.precision(15);
	
	fout_p << D << endl;

	for (unsigned int n = 0; n < D; ++n) {
			
		fout_p << sorted_part_poles(n) << endl;

		for(unsigned int itr_comp = 0; itr_comp < D; itr_comp++) {

			fout_p << sorted_part_residues(n, itr_comp) << "  ";
		}

		fout_p << endl;
	}


	fout_p.clear();
	fout_p.close();

	//cout << "The OpRS propagator has been printed!" << endl;
}





