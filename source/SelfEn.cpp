#include <omp.h>

#include "SelfEn.h"


using namespace std;


M_II::M_II(const SpProp& G_ref): D( G_ref.X_np.cols() ), n_max( G_ref.X_np.rows() ), k_max( G_ref.Y_kp.rows() ), M( n_max*n_max*k_max, D ) {

	M.setZero();
}

M_II::~M_II() {

}






N_II::N_II(const SpProp& G_ref): D( G_ref.X_np.cols() ), n_max( G_ref.X_np.rows() ), k_max( G_ref.Y_kp.rows() ), N( D, k_max*k_max*n_max ) {

	N.setZero();
}

N_II::~N_II() {

}







E_gr_II::E_gr_II(const SpProp& G_ref): n_max( G_ref.X_np.rows() ), k_max( G_ref.Y_kp.rows() ), E( n_max*n_max*k_max, n_max*n_max*k_max ) {

	E.setZero();
}

E_gr_II::~E_gr_II() {

}







E_less_II::E_less_II(const SpProp& G_ref): n_max( G_ref.X_np.rows() ), k_max( G_ref.Y_kp.rows() ), E( k_max*k_max*n_max, k_max*k_max*n_max ) {

	E.setZero();
}

E_less_II::~E_less_II() {

}








SelfEn::SelfEn(const SpBasis& SpBas, const SpProp& G, const double pair): Sigma_inf(SpBas.Num_Levels), M_2p1h(G), N_2h1p(G), E_2p1h(G), E_2h1p(G), g(pair) {

	Sigma_inf.setZero();
}

	
SelfEn::~SelfEn() {

}


void SelfEn::BuildStatic(const SpProp& G) {
	
	//cout << "Building static self-energy..." << endl;

	Sigma_inf.setZero();

	for (int itr_alpha = 0; itr_alpha < G.Y_kp.cols(); ++itr_alpha) {
	
		for (int itr_h = 0; itr_h < G.Y_kp.rows(); ++itr_h) {

			Sigma_inf(itr_alpha) += pow(G.Y_kp(itr_h, itr_alpha) , 2);
		}

		Sigma_inf(itr_alpha) *= -g/2.;
	}

	//cout << "The static self-energy has been built!" << endl;
}


void SelfEn::Tadpole(const SpBasis& SpBas, const string path, const string id) {

	int num_pairs = SpBas.Num_Holes/2;

 	DenseVec poles = SpBas.e_0 + Sigma_inf;
   	DenseMat residues(SpBas.Num_Levels, SpBas.Num_Levels);

	for (int itr_alpha = 0; itr_alpha < SpBas.Num_Levels; ++itr_alpha) {

		residues(itr_alpha, itr_alpha) = 1.;
	}


	//cout << "The static self-energy has been diagonalized! Now printing it..." << endl;

	ofstream fout_h(path + "spprop_backward_" + id + ".dat");
	fout_h.precision(15);
	fout_h << num_pairs << endl;

	ofstream fout_p(path + "spprop_forward_" + id + ".dat");
	fout_p.precision(15);
	fout_p << SpBas.Num_Levels - num_pairs << endl;


	for (int itr_pole = 0; itr_pole < residues.rows(); ++itr_pole) {

		if (itr_pole < num_pairs) {

			fout_h << poles(itr_pole) << endl;

			for (int itr_alpha = 0; itr_alpha < residues.cols(); ++itr_alpha) {

				fout_h << residues(itr_pole, itr_alpha) << "  ";
			}
		
			fout_h << endl;
		}
		else {

			fout_p << poles(itr_pole) << endl;

			for (int itr_alpha = 0; itr_alpha < residues.cols(); ++itr_alpha) {

				fout_p << residues(itr_pole, itr_alpha) << "  ";
			}
		
			fout_p << endl;
		}
	}

	//cout << "The Static self-energy has been printed!" << endl;
}




void SelfEn::PrintDyson(const SpBasis& SpBas, DenseVec& eigenvalues, DenseMat& eigenvectors, const string path, const string id) {
	
	unsigned int num_pairs = SpBas.Num_Holes/2;
	unsigned int num_rows = eigenvectors.rows();
	unsigned int num_cols = eigenvectors.cols();

	
	DenseMat Z_i(num_rows, num_cols);
	DenseVec esp_i(num_rows);


	unsigned int row = 0;

	for (unsigned int itr_row = 0; itr_row < num_rows; ++itr_row) {

		Z_i.row(row) = eigenvectors.row(itr_row);
		esp_i(row) = eigenvalues(itr_row);
		++row;
	}

	//cout << "Printing new propagator..." << endl;

	ofstream fout_h(path + "spprop_backward_" + id + ".dat");
	fout_h.precision(15);

	fout_h << "1               " << endl; 


	ofstream fout_p(path + "spprop_forward_" + id + ".dat");
	fout_p.precision(15);
	
	fout_p << "1               " << endl; 


	double num_filled = 0.;
	double num_filled_prev = 0.;


	int num_h = 0;

	for(unsigned int itr_esp_i = 0; itr_esp_i < num_rows; ++itr_esp_i) {

		num_filled_prev = num_filled;

		for(unsigned int itr_comp = 0; itr_comp < num_cols; ++itr_comp) {

			num_filled += pow( Z_i(itr_esp_i, itr_comp), 2 );
		}

		if( abs(num_filled - num_pairs) < abs(num_filled_prev - num_pairs) ) {
			
			fout_h << esp_i(itr_esp_i) << endl;

			for(unsigned int itr_comp = 0; itr_comp < num_cols; ++itr_comp) {

				fout_h << Z_i(itr_esp_i, itr_comp) << "  ";
			}

			fout_h << endl;
			
			++num_h;
		}
		else {

			fout_p << esp_i(itr_esp_i) << endl;

			for(unsigned int itr_comp = 0; itr_comp < num_cols; itr_comp++) {

				fout_p << Z_i(itr_esp_i, itr_comp) << "  ";
			}

			fout_p << endl;
		}
	}

	int num_p = num_rows - num_h;

	fout_h.seekp(0);
	fout_h << num_h << endl;

	fout_h.clear();
	fout_h.close();


	fout_p.seekp(0);
	fout_p << num_p << endl;;

	fout_p.clear();
	fout_p.close();

	//cout << "The propagator has been printed!" << endl;
}



void SelfEn::ChargeDysMat(const SpBasis& SpBas, const unsigned int alpha) {

	DenseMat DenseDysMat( 1 + M_2p1h.M.rows() + N_2h1p.N.cols(), 1 + M_2p1h.M.rows() + N_2h1p.N.cols() );

	DenseDysMat << Sigma_inf(alpha) + SpBas.e_0(alpha), M_2p1h.M.col(alpha).transpose(), N_2h1p.N.row(alpha),
					M_2p1h.M.col(alpha), E_2p1h.E, DenseMat::Zero(M_2p1h.M.rows(), N_2h1p.N.cols()),
					N_2h1p.N.row(alpha).transpose(), DenseMat::Zero( N_2h1p.N.cols(), M_2p1h.M.rows()), E_2h1p.E;

	//cout << DenseDysMat << endl;

	// Prune the matrix from diagonal rows and columns (if they aren't the first row/column)

	vector <unsigned int> to_keep;

	to_keep.push_back(0);

	for (int itr_diag = 1; itr_diag < DenseDysMat.rows(); ++itr_diag) {

		bool nonzero = false;

		for (int itr_row = 0; itr_row < DenseDysMat.rows(); ++itr_row) {

			if ( abs(DenseDysMat(itr_row, itr_diag)) > 10E-7 && itr_row != itr_diag) {

				to_keep.push_back(itr_diag);
				nonzero = true;
				break;
			}
		}
	}


	DenseMat DenseDysMat_reduced(to_keep.size(), to_keep.size());
	DenseDysMat_reduced.setZero();


	for (int itr_row = 0; itr_row < to_keep.size(); ++itr_row) {

		for (int itr_col = 0; itr_col < to_keep.size(); ++itr_col) {

			DenseDysMat_reduced(itr_row, itr_col) = DenseDysMat(to_keep[itr_row], to_keep[itr_col]);
		}
	}

	//cout << DenseDysMat_reduced << endl;
	//cout << endl;
	//cout << endl;

	//cout << "alpha = " << alpha << endl;
	//cout << M_2p1h.M.col(alpha) << endl;
	//cout << endl;

	//cout << "The Dyson matrix has dimension " << DenseDysMat_reduced.rows() << "x" << DenseDysMat_reduced.cols() << endl;
	
	DysMat = DenseDysMat_reduced.sparseView().pruned(10E-7);
}




void SelfEn::DiagonalizeDysMat(const SpBasis& SpBas, const SpProp& G, const string path, const string id) {

	DenseVec Poles(0);
	DenseMat Residues(0, SpBas.Num_Levels);

	for (int alpha = 0; alpha < SpBas.Num_Levels; ++alpha) {

		ChargeDysMat(SpBas, alpha);

		Eigen::SelfAdjointEigenSolver < Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > solver(DysMat);

 		DenseVec eigenvalues = solver.eigenvalues();
   		DenseMat eigenvectors = solver.eigenvectors();

		Poles.conservativeResize(Poles.size() + eigenvalues.size());
		Poles.tail(eigenvalues.size()) = eigenvalues;

		int old_num_rows = Residues.rows();

		Residues.conservativeResizeLike( Eigen::MatrixXd::Zero( Residues.rows() + eigenvalues.size(), Residues.cols()) );

		int num_non_zero = 0;

		for (int itr_new_poles = 0; itr_new_poles < eigenvalues.size(); ++itr_new_poles) {

			if ( abs(eigenvectors(0, itr_new_poles)) > 10E-7 ) {

				Residues( old_num_rows + num_non_zero, alpha) = eigenvectors(0, itr_new_poles);
				Poles(old_num_rows + num_non_zero) = eigenvalues(itr_new_poles);
				++num_non_zero;
			}
		}

		Residues.conservativeResize( old_num_rows + num_non_zero, Residues.cols() );
		Poles.conservativeResize( old_num_rows + num_non_zero );

		//cout << "The self-energy has been diagonalized for " << "alpha = " << alpha << "..." << endl;
	}

	// Ordering the Poles and Residues 

	DenseVec sorted_Poles(Poles);
	DenseMat sorted_Residues(Residues.rows(), SpBas.Num_Levels);

	sort(sorted_Poles.data(), sorted_Poles.data() + sorted_Poles.size());

	for (int itr = 0; itr < sorted_Poles.size(); ++itr) {
        
		Eigen::Index index;
		(Poles.array() - sorted_Poles(itr)).abs().minCoeff(&index); // Find index of closest element
        sorted_Residues.row(itr) = Residues.row(index); // Rearrange rows
    }

	PrintDyson(SpBas, sorted_Poles, sorted_Residues, path, id);
}



void SelfEn::sc0(const SpBasis& SpBas, SpProp& G) {

	string path = "./SpProp/";

	double E0_GS(0.);
	double E1_GS(1.);
	double E_ref(G.MGK(SpBas));
	
	unsigned int limit = 30;
	unsigned int count = 0;

	while ( abs(E1_GS - E0_GS) > 10E-8 && count < limit) {

		//cout << "Entering sc0 iterations..." << endl;
		//cout << "Ground state energy = " << E_ref << "     Difference wrt to reference energy = " << 0. << endl;

		string id = "sc0_itr" + to_string(count);
		
		DiagonalizeDysMat(SpBas, G, path, id);

		G.ChargeNewProp(SpBas, path, id);
		BuildStatic(G);

		E0_GS = E1_GS;
		E1_GS = G.MGK(SpBas);

		//cout << endl;
		++count;
		
		//cout << "Iter " << count << ": Ground state energy = " << E1_GS << "     Difference wrt to reference energy = " << E1_GS - E_ref << endl;
	}

	if (count-- == limit ) {

		cout << "The sc0 iterations exited before convergence...   Diff_En = " << abs(E1_GS - E0_GS) << endl; 
	}

	string oldname = path + "spprop_backward_sc0_itr" + to_string(count) + ".dat";
	string newname = path + "spprop_backward_sc0.dat";
		
	rename(oldname.c_str(), newname.c_str());

	oldname = path + "spprop_forward_sc0_itr" + to_string(count) + ".dat";
	newname = path + "spprop_forward_sc0.dat";

	rename(oldname.c_str(), newname.c_str());
}



void SelfEn::BuildADC2(const SpProp& G) {

	M_2p1h.M.setZero();
	E_2p1h.E.setZero();

	// Computing M_2p1h and E_2p1h

	for (int n1 = 0; n1 < G.X_np.rows(); ++n1) {
		
		for (int n2 = 0; n2 < G.X_np.rows(); ++n2) {
	
			for (int k3 = 0; k3 < G.Y_kp.rows(); ++k3) {

				for (int alpha = 0; alpha < G.X_np.cols(); ++alpha) {
					
					for (int gamma = 0; gamma < G.X_np.cols(); ++gamma) {

						M_2p1h(n1, n2, k3, alpha) += G.X_np(n1, gamma)*G.X_np(n2, gamma)*G.Y_kp(k3, alpha); 
					}

					M_2p1h(n1, n2, k3, alpha) *= -g/2.;
				}

				E_2p1h(n1, n2, k3, n1, n2, k3) = G.esp_n(n1) + G.esp_n(n2) - G.esp_k(k3);
			}
		}
	}



	N_2h1p.N.setZero();
	E_2h1p.E.setZero();


	// Computing N_2h1p and E_2h1p

	for (int k1 = 0; k1 < G.Y_kp.rows(); ++k1) {
		
		for (int k2 = 0; k2 < G.Y_kp.rows(); ++k2) {
	
			for (int n3 = 0; n3 < G.X_np.rows(); ++n3) {

				for (int alpha = 0; alpha < G.X_np.cols(); ++alpha) {
				
					for (int gamma = 0; gamma < G.Y_kp.cols(); ++gamma) {

						N_2h1p(alpha, k1, k2, n3) += G.Y_kp(k1, gamma)*G.Y_kp(k2, gamma)*G.X_np(n3, alpha); 
					}

					N_2h1p(alpha, k1, k2, n3) *= -g/2.;
				}

				E_2h1p(k1, k2, n3, k1, k2, n3) = G.esp_k(k1) + G.esp_k(k2) - G.esp_n(n3);
			}
		}
	}
}




void SelfEn::BuildTDA(const SpProp& G) {

	// Computing M_2p1h and E_2p1h
	
	M_2p1h.M.setZero();
	E_2p1h.E.setZero();

	for (int n1 = 0; n1 < G.X_np.rows(); ++n1) {
		
		for (int n2 = 0; n2 < G.X_np.rows(); ++n2) {
	
			for (int k3 = 0; k3 < G.Y_kp.rows(); ++k3) {

				E_2p1h(n1, n2, k3, n1, n2, k3) = G.esp_n(n1) + G.esp_n(n2) - G.esp_k(k3);

				for (int alpha = 0; alpha < G.X_np.cols(); ++alpha) {
					
					for (int gamma = 0; gamma < G.X_np.cols(); ++gamma) {

						M_2p1h(n1, n2, k3, alpha) += G.X_np(n1, gamma)*G.X_np(n2, gamma)*G.Y_kp(k3, alpha); 

						for (int n4 = 0; n4 < G.X_np.rows(); ++n4) {
		
							for (int n5 = 0; n5 < G.X_np.rows(); ++n5) {
				
								E_2p1h(n1, n2, k3, n4, n5, k3) += -g/2.*G.X_np(n1, alpha)*G.X_np(n2, alpha)*G.X_np(n4, gamma)*G.X_np(n5, gamma);
							}
						}
					}
					
					M_2p1h(n1, n2, k3, alpha) *= -g/2.;
				}
			}
		}
	}



	// Computing N_2h1p and E_2h1p

	N_2h1p.N.setZero();
	E_2h1p.E.setZero();

	for (int k1 = 0; k1 < G.Y_kp.rows(); ++k1) {
		
		for (int k2 = 0; k2 < G.Y_kp.rows(); ++k2) {
	
			for (int n3 = 0; n3 < G.X_np.rows(); ++n3) {

				E_2h1p(k1, k2, n3, k1, k2, n3) = G.esp_k(k1) + G.esp_k(k2) - G.esp_n(n3);
				
				for (int alpha = 0; alpha < G.X_np.cols(); ++alpha) {
				
					for (int gamma = 0; gamma < G.Y_kp.cols(); ++gamma) {

						N_2h1p(alpha, k1, k2, n3) += G.Y_kp(k1, gamma)*G.Y_kp(k2, gamma)*G.X_np(n3, alpha); 
						
						for (int k4 = 0; k4 < G.Y_kp.rows(); ++k4) {
		
							for (int k5 = 0; k5 < G.Y_kp.rows(); ++k5) {
				
								E_2h1p(k1, k2, n3, k4, k5, n3) += g/2.*G.Y_kp(k1, alpha)*G.Y_kp(k2, alpha)*G.Y_kp(k4, gamma)*G.Y_kp(k5, gamma);
							}
						}
					}

					N_2h1p(alpha, k1, k2, n3) *= -g/2.;
				}
			}
		}
	}
}



/* Deprecated*/
/*
void SelfEn::BuildADC3(const SpProp& G) {

	BuildTDA(G);

	for (int n1 = 0; n1 < G.X_np.rows(); ++n1) {
		
		for (int n2 = 0; n2 < G.X_np.rows(); ++n2) {
	
			for (int k3 = 0; k3 < G.Y_kp.rows(); ++k3) {

				for (int alpha = 0; alpha < G.Y_kp.cols(); ++alpha) {

					for (int k4 = 0; k4 < G.Y_kp.rows(); ++k4) {

						for (int k5 = 0; k5 < G.Y_kp.rows(); ++k5) {

							for (int gamma = 0; gamma < G.Y_kp.cols(); ++gamma) {

								for (int delta = 0; delta < G.Y_kp.cols(); ++delta) {

									M_2p1h(n1, n2, k3, alpha) += pow(g/2., 2)*pow(G.Y_kp(k4, gamma), 2)*pow(G.Y_kp(k5, gamma), 2)*G.X_np(n1, delta)*G.X_np(n2, delta)*G.Y_kp(k3, alpha)/(G.esp_k(k4) + G.esp_k(k5) - G.esp_n(n1) - G.esp_n(n2));
								}
							}
						}
					}
				}
			}
		}
	}


	for (int k1 = 0; k1 < G.Y_kp.rows(); ++k1) {
		
		for (int k2 = 0; k2 < G.Y_kp.rows(); ++k2) {
	
			for (int n3 = 0; n3 < G.X_np.rows(); ++n3) {

				for (int alpha = 0; alpha < G.Y_kp.cols(); ++alpha) {

					for (int n4 = 0; n4 < G.X_np.rows(); ++n4) {

						for (int n5 = 0; n5 < G.X_np.rows(); ++n5) {

							for (int gamma = 0; gamma < G.Y_kp.cols(); ++gamma) {

								for (int delta = 0; delta < G.Y_kp.cols(); ++delta) {

									N_2h1p(alpha, k1, k2, n3) += pow(g/2., 2)*pow(G.X_np(n4, gamma), 2)*pow(G.X_np(n5, gamma), 2)*G.Y_kp(k1, delta)*G.Y_kp(k2, delta)*G.X_np(n3, alpha)/(G.esp_k(k1) + G.esp_k(k2) - G.esp_n(n4) - G.esp_n(n5));
								}
							}
						}
					}
				}
			}
		}
	}
}
*/



void SelfEn::BuildADC3(const SpProp& G) {

	BuildTDA(G);

	// Build t0 matrix

	Eigen::Tensor<double, 4> t(G.X_np.rows(), G.X_np.rows(), G.Y_kp.rows(), G.Y_kp.rows());
	t.setZero();

	for (int n1 = 0; n1 < G.X_np.rows(); ++n1) {
		
		for (int n2 = 0; n2 < G.X_np.rows(); ++n2) {
	
			for (int k3 = 0; k3 < G.Y_kp.rows(); ++k3) {

				for (int k4 = 0; k4 < G.Y_kp.rows(); ++k4) {

					double one_over_e = 1./(G.esp_k(k3) + G.esp_k(k4) - G.esp_n(n1) - G.esp_n(n2));

					for (int alpha = 0; alpha < G.Y_kp.cols(); ++alpha) {

						for (int beta = 0; beta < G.Y_kp.cols(); ++beta) {
						
							t(n1, n2, k3, k4) += G.X_np(n1, alpha)*G.X_np(n2, alpha)*G.Y_kp(k3, beta)*G.Y_kp(k4, beta);
						}
					}

					t(n1, n2, k3, k4) *= -g/2.*one_over_e;
				}
			}
		}
	}

	///// Now calculate M_II and N_II

	for (int n1 = 0; n1 < G.X_np.rows(); ++n1) {
		
		for (int n2 = 0; n2 < G.X_np.rows(); ++n2) {
	
			for (int k3 = 0; k3 < G.Y_kp.rows(); ++k3) {

				for (int alpha = 0; alpha < G.Y_kp.cols(); ++alpha) {

					for (int k4 = 0; k4 < G.Y_kp.rows(); ++k4) {

						for (int k5 = 0; k5 < G.Y_kp.rows(); ++k5) {

							for (int gamma = 0; gamma < G.Y_kp.cols(); ++gamma) {

								M_2p1h(n1, n2, k3, alpha) += -g/2*t(n1, n2, k4, k5)*G.Y_kp(k4, gamma)*G.Y_kp(k5, gamma)*G.Y_kp(k3, alpha);
							}
						}
					}
				}
			}
		}
	}

	for (int k1 = 0; k1 < G.Y_kp.rows(); ++k1) {
		
		for (int k2 = 0; k2 < G.Y_kp.rows(); ++k2) {
	
			for (int n3 = 0; n3 < G.X_np.rows(); ++n3) {

				for (int alpha = 0; alpha < G.Y_kp.cols(); ++alpha) {

					for (int n4 = 0; n4 < G.X_np.rows(); ++n4) {

						for (int n5 = 0; n5 < G.X_np.rows(); ++n5) {

							for (int gamma = 0; gamma < G.Y_kp.cols(); ++gamma) {

								N_2h1p(alpha, k1, k2, n3) += -g/2.*t(n4, n5, k1, k2)*G.X_np(n4, gamma)*G.X_np(n5, gamma)*G.X_np(n3, alpha);	
							}
						}
					}
				}
			}
		}
	}
}



Eigen::Tensor<double, 8> SelfEn:: IterateCCD(const Eigen::Tensor<double, 8>& t_0, const double g_proxy, const SpProp& G, const unsigned int P, const double alpha) {

	const unsigned int D = G.X_np.cols();
	unsigned int n_max = (D - P);
	unsigned int k_max = P;
	
	
	Eigen::Tensor<double, 4> one_over_esp_ab_ij(n_max, n_max, k_max, k_max);
	one_over_esp_ab_ij.setZero();

	for (unsigned int a = 0; a < n_max; ++a) {

		for (unsigned int b = 0; b < n_max; ++b) {
		
			for (unsigned int i = 0; i < k_max; ++i) {

				for (unsigned int j = 0; j < k_max; ++j) {

					one_over_esp_ab_ij(a, b, i, j) = 1./( G.esp_k(i) + G.esp_k(j) - G.esp_n(a) - G.esp_n(b));
				}
			}
		}
	}

	auto v = [g_proxy](unsigned int s_1, unsigned int s_2, unsigned int s_3, unsigned int s_4) { 
		
		if (s_1 != s_2 && s_1 == s_3 && s_2 == s_4) {

			return -g_proxy/2;
		}
		else if (s_1 != s_2 && s_1 != s_3 && s_2 != s_4) {

			return g_proxy/2.;
		}
		
		return  0.; 
	};

	double tol = 1E-7;

	auto NormDiff = [](const auto& A, const auto& B) {
    
		auto T = A - B;
		const Eigen::Tensor<double, 0> frob_norm_tens = T.square().sum().sqrt();
  		const double frob_norm = frob_norm_tens.coeff();

		return frob_norm;
	};

	Eigen::Tensor<double, 8> t(t_0);
	Eigen::Tensor<double, 8> t_n(n_max, 2, n_max, 2, k_max, 2, k_max, 2);

	while ( NormDiff(t, t_n) > tol) {

		t_n.setZero();

		// First term
		for (unsigned int a = 0; a < n_max; ++a) {

			for (unsigned int i = 0; i < k_max; ++i) {

				t_n(a, 0, a, 1, i, 0, i, 1) += -g_proxy/2;
			}
		}

		
		// Second term
		for (unsigned int a = 0; a < n_max; ++a) {

			for (unsigned int i = 0; i < k_max; ++i) {

				for (unsigned int j = 0; j < k_max; ++j) {

					for (unsigned int c = 0; c < n_max; ++c) {

						for (unsigned int s_c = 0; s_c < 2; ++s_c) {
							
							for (unsigned int s_d = 0; s_d < 2; ++s_d) {

								t_n(a, 0, a, 1, i, 0, j, 1) += 0.5*v(0, 1, s_c, s_d)*t(c, s_c, c, s_d, i, 0, j, 1);
							}
						}
					}
				}
			}
		}
		

		// Third term
		for (unsigned int a = 0; a < n_max; ++a) {

			for (unsigned int b = 0; b < n_max; ++b) {

				for (unsigned int i = 0; i < k_max; ++i) {

					for (unsigned int k = 0; k < k_max; ++k) {

						for (unsigned int s_k = 0; s_k < 2; ++s_k) {
							
							for (unsigned int s_l = 0; s_l < 2; ++s_l) {

								t_n(a, 0, b, 1, i, 0, i, 1) += .5*v(s_k, s_l, 0, 1)*t(a, 0, b, 1, k, s_k, k, s_l);
							}
						}
					}
				}
			}
		}

		// Other terms
		for (unsigned int a = 0; a < n_max; ++a) {

			for (unsigned int b = 0; b < n_max; ++b) {

				for (unsigned int i = 0; i < k_max; ++i) {

					for (unsigned int j = 0; j < k_max; ++j) {

						for (unsigned int k = 0; k < k_max; ++k) {

							for (unsigned int s_k = 0; s_k < 2; ++s_k) {
							
								for (unsigned int s_l = 0; s_l < 2; ++s_l) {

									for (unsigned int c = 0; c < n_max; ++c) {

										for (unsigned int s_c = 0; s_c < 2; ++s_c) {
							
											for (unsigned int s_d = 0; s_d < 2; ++s_d) {

												// Fourth term
												t_n(a, 0, b, 1, i, 0, j, 1) += .25*v(s_k, s_l, s_c, s_d)*t(c, s_c, c, s_d, i, 0, j, 1)*t(a, 0, b, 1, k, s_k, k, s_l);
												
												// Fifth term
												t_n(a, 0, b, 1, i, 0, j, 1) += v(s_k, s_l, s_c, s_d)*t(a, 0, c, s_c, i, 0, k, s_k)*t(b, 1, c, s_d, j, 1, k, s_l);
												t_n(a, 0, b, 1, i, 0, j, 1) += -v(s_k, s_l, s_c, s_d)*t(a, 0, c, s_c, j, 1, k, s_k)*t(b, 1, c, s_d, i, 0, k, s_l);

												// Sixth term
												t_n(a, 0, b, 1, i, 0, j, 1) += -.5*v(s_k, s_l, s_c, s_d)*t(c, s_d, c, s_c, i, 0, k, s_k)*t(a, 0, b, 1, k, s_l, j, 1);
												t_n(a, 0, b, 1, i, 0, j, 1) += .5*v(s_k, s_l, s_c, s_d)*t(c, s_d, c, s_c, j, 1, k, s_k)*t(a, 0, b, 1, k, s_l, i, 0);

												// Seventh term
												t_n(a, 0, b, 1, i, 0, j, 1) += -.5*v(s_k, s_l, s_c, s_d)*t(a, 0, c, s_c, k, s_l, k, s_k)*t(c, s_d, b, 1, i, 0, j, 1);
												t_n(a, 0, b, 1, i, 0, j, 1) += .5*v(s_k, s_l, s_c, s_d)*t(b, 1, c, s_c, k, s_l, k, s_k)*t(c, s_d, a, 0, i, 0, j, 1);;

											}
										}
									}
								}
							}
						}
					}
				}
			}
		}


		// Prefactor
		for (unsigned int a = 0; a < n_max; ++a) {

			for (unsigned int s_a = 0; s_a < 2; ++s_a) {

				for (unsigned int b = 0; b < n_max; ++b) {

					for (unsigned int s_b = s_a + 1; s_b < 2; ++s_b) {

						for (unsigned int i = 0; i < k_max; ++i) {

							for (unsigned int s_i = 0; s_i < 2; ++s_i) {

								for (unsigned int j = 0; j < k_max; ++j) {
							
									for (unsigned int s_j = s_i + 1; s_j < 2; ++s_j) {

										t_n(a, s_a, b, s_b, i, s_i, j, s_j) *= one_over_esp_ab_ij(a, b, i, j);
									}
								}
							}
						}
					}
				}
			}
		}


		// Fix antisymmetry
		for (unsigned int a = 0; a < n_max; ++a) {

			for (unsigned int b = 0; b < n_max; ++b) {

				for (unsigned int i = 0; i < k_max; ++i) {

					for (unsigned int j = 0; j < k_max; ++j) {

						t_n(a, 1, b, 0, i, 0, j, 1) = -t_n(b, 0, a, 1, i, 0, j, 1);
						t_n(a, 0, b, 1, i, 1, j, 0) = -t_n(a, 0, b, 1, j, 0, i, 1);
						t_n(a, 1, b, 0, i, 1, j, 0) = t_n(b, 0, a, 1, j, 0, i, 1);
					}
				}
			}
		}

		t = alpha*t_n + (1. - alpha)*t;

		if (NormDiff(t, t_n) > 1.E20) {

			cerr << "The calculation of t is diverging for g = " << g_proxy << endl;
			exit(5);
		}

		//cout << NormDiff(t, t_n) << endl;
	}

	cout << "t_ij^ab amplitude converged to " << tol << " for g = " << g_proxy << endl;
	
	return t;
}




void SelfEn::BuildADC3D(const SpProp& G, const unsigned int P) {

	BuildTDA(G);

	const unsigned int D = G.X_np.cols();
	
	// Build t0 matrix

	unsigned int n_max = (D - P);
	unsigned int k_max = P;

	Eigen::Tensor<double, 8> t_0(n_max, 2, n_max, 2, k_max, 2, k_max, 2);
	t_0.setZero();
	
	Eigen::Tensor<double, 2> one_over_esp_ab_ij(n_max, k_max);
	one_over_esp_ab_ij.setZero();

	for (unsigned int a = 0; a < n_max; ++a) {

		for (unsigned int i = 0; i < k_max; ++i) {

			one_over_esp_ab_ij(a, i) = 1./( G.esp_k(i) + G.esp_k(i) - G.esp_n(a) - G.esp_n(a) );
		}
	}

	for (unsigned int a = 0; a < n_max; ++a) {
		
		for (unsigned int i = 0; i < k_max; ++i) {

			double g_one_over_e = -g/2.*one_over_esp_ab_ij(a, i);
	
			for (unsigned int s = 0; s < 2; ++s) {

				t_0(a, s, a, 1-s, i, s, i, 1-s) = g_one_over_e;
				t_0(a, s, a, 1-s, i, 1-s, i, s) += -g_one_over_e;
			}
		}
	}

	// t_0 built

	//double g_proxy = 1.273;
	//auto t_1 = IterateCCD(t_0, g_proxy, G, P, 0.4);
	auto t = IterateCCD(t_0, g, G, P, 0.4);


	///// Now calculate M_II and N_II

	for (int alpha = 0; alpha < P; ++alpha) {

		for (int n1 = 0; n1 < n_max; ++n1) {

			for (int n2 = 0; n2 < n_max; ++n2) {

				for (int k4 = 0; k4 < k_max; ++k4) {

					M_2p1h(n1, n2, alpha, alpha) += -g/2.*t(n1, 1, n2, 0, k4, 1, k4, 0);
				}
			}
		}
	}

	for (int alpha = P; alpha < D; ++alpha) {
	
		for (int k1 = 0; k1 < k_max; ++k1) {

			for (int k2 = 0; k2 < k_max; ++k2) {

				for (int n4 = 0; n4 < G.X_np.rows(); ++n4) {

					N_2h1p(alpha, k1, k2, alpha - P) += -g/2.*t(n4, 1, n4, 0, k1, 1, k2, 0);	
				}
			}
		}
	}
}