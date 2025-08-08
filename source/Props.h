#ifndef __Props_h__
#define __Props_h__

#include <cmath>
#include <cstring>

#include "SelfEn.h"
#include "Basis.h"

using namespace std;


struct EdgePropr;


class SpProp {

	public:
	
	SpProp();
	~SpProp();
	SpProp(const SpProp&);
	SpProp& operator = (const SpProp&);
	
	void Build_G0(const SpBasis&);
	
	void Print();
	void ChargeNewProp(const SpBasis&, const string, const string);
	void PrintFile(const SpBasis&, const string, const string);
	void SpecFunc(const string, const double);
	double MGK(const SpBasis&) const;
	void GenerateOpRS();
	
	
	public:
	
	DenseMat Y_kp;  
	DenseMat X_np;

	DenseVec esp_k;
	DenseVec esp_n;
};

#endif
