#ifndef EIGENVALUE_H 
#define EIGENVALUE_H 

#include "tnt/jama_eig.h"

//////////////////////////////////////////////////////////////////////
// eigenvalues of 3x3 non-symmetric matrix
//////////////////////////////////////////////////////////////////////
void inline computeEigenvalues3x3(
		float dout[3], 
		float a[3][3])
{
	TNT::Array2D<float> A = TNT::Array2D<float>(3,3, &a[0][0]);
	TNT::Array1D<float> eig = TNT::Array1D<float>(3);
	TNT::Array1D<float> eigImag = TNT::Array1D<float>(3);
	JAMA::Eigenvalue<float> jeig = JAMA::Eigenvalue<float>(A);
	jeig.getRealEigenvalues(eig);

	// complex ones
	jeig.getImagEigenvalues(eigImag);
	dout[0]  = sqrt(eig[0]*eig[0] + eigImag[0]*eigImag[0]);
	dout[1]  = sqrt(eig[1]*eig[1] + eigImag[1]*eigImag[1]);
	dout[2]  = sqrt(eig[2]*eig[2] + eigImag[2]*eigImag[2]);
}

#undef rfabs 
#undef ROT
#endif
