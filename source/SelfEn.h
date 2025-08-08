#ifndef __SelfEn_h__
#define __SelfEn_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cmath>

#include "Basis.h"
#include "Props.h"
#include "Functions.h"



class SpProp;

class M_II {

	public:

	M_II(const SpProp&);
	~M_II();

	double& operator() (const unsigned int n1, const unsigned int n2, const unsigned int k3, const unsigned int alpha) {

		return M(n1*n_max*k_max + n2*k_max + k3, alpha);
	};


	double& operator() (const unsigned int n, const unsigned int k, const unsigned int alpha) {

		return M(n*k_max + k, alpha);
	};


	public:

	unsigned int D;
	unsigned int n_max;
	unsigned int k_max;

	DenseMat M;
};



class N_II {

	public:

	N_II(const SpProp&);
	~N_II();


	double& operator() (const unsigned int alpha, const unsigned int k1, const unsigned int k2, const unsigned int n3) {

		return N(alpha, k1*k_max*n_max + k2*n_max + n3);
	};

	double& operator() (const unsigned int alpha, const unsigned int k, const unsigned int n) {

		return N(alpha, k*n_max + n);
	};

	public:

	unsigned int D;
	unsigned int n_max;
	unsigned int k_max;

	DenseMat N;
};



class E_gr_II {

	public:

	E_gr_II(const SpProp&);
	~E_gr_II();

	double& operator() (const unsigned int n1, const unsigned int n2, const unsigned int k3, const unsigned int n4, const unsigned int n5, const unsigned int k6) {

		return E(n1*n_max*k_max + n2*k_max + k3, n4*n_max*k_max + n5*k_max + k6);
	};

	double& operator() (const unsigned int n, const unsigned int k, const unsigned int n_p, const unsigned int k_p) {

		return E(n*k_max + k, n_p*k_max + k_p);
	};
	

	public:

	unsigned int n_max;
	unsigned int k_max;

	DenseMat E;
};


class E_less_II {

	public:

	E_less_II(const SpProp&);
	~E_less_II();

	double& operator() (const unsigned int k1, const unsigned int k2, const unsigned int n3, const unsigned int k4, const unsigned int k5, const unsigned int n6) {

		return E(k1*k_max*n_max + k2*n_max + n3, k4*k_max*n_max + k5*n_max + n6);
	};

	double& operator() (const unsigned int k, const unsigned int n, const unsigned int k_p, const unsigned int n_p) {

		return E(k*n_max + n, k_p*n_max + n_p);
	};


	public:

	unsigned int n_max;
	unsigned int k_max;

	DenseMat E;
};






class SelfEn {

	public:
	
	SelfEn(const SpBasis&, const SpProp& G, const double); 
	~SelfEn();
	void BuildStatic(const SpProp&);
	void Tadpole(const SpBasis&, const string, const string);
	void PrintDyson(const SpBasis&, DenseVec&, DenseMat&, const string, const string);
	void ChargeDysMat(const SpBasis&, const unsigned int);
	void DiagonalizeDysMat(const SpBasis&, const SpProp&, const string, const string);
	void sc0(const SpBasis&, SpProp&);
	void BuildADC2(const SpProp&);
	void BuildTDA(const SpProp&);
	void BuildADC3(const SpProp&);
	Eigen::Tensor<double, 8> IterateCCD(const Eigen::Tensor<double, 8>&, const double, const SpProp&, const unsigned int);
	void BuildADC3D(const SpProp&, const unsigned int);

	
	public:

	const double g;

	DenseVec Sigma_inf;
	SparseMat DysMat;

	DenseVec E_gr;
	DenseVec E_less;

	M_II M_2p1h;
	N_II N_2h1p;

	E_gr_II E_2p1h;
	E_less_II E_2h1p;
};

#endif
