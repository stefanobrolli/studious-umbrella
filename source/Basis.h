#ifndef __Basis_h
#define __Basis_h

#include <complex>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <eigen3/Eigen/Eigenvalues>


using namespace std;

const complex <double> i(0, 1.);

typedef Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> DenseMat;
typedef Eigen::SparseMatrix<double> SparseMat;
typedef Eigen::Matrix <complex <double>, Eigen::Dynamic, Eigen::Dynamic> DenseCompMat;
typedef Eigen::Vector <double, Eigen::Dynamic> DenseVec;
typedef Eigen::Vector < complex <double>, Eigen::Dynamic> DenseCompVec;


class SpBasis {

	public:
	
	SpBasis(const unsigned int, const unsigned int, double);
	~SpBasis();

	
	public:
	
	unsigned int Num_Levels;
	unsigned int Num_Holes;
	DenseVec e_0;
};

#endif
