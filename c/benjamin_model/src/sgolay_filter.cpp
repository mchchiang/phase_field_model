// sgolay_filter.cpp

#include <iostream>
#include <cmath>
#include <vector>
#include <armadillo>
#include "filter.h"

using std::vector;
using namespace arma;

class SgolayFilter {

private:
  int degree;
  int npoints;
  vector<double> coeffs;
  vector<double> derivCoeffs;

public:
  SgolayFilter(int z, int npts) {
    degree = z;
    npoints = npts;
    
    // Find the coefficients
    coeffs = vector<double>(npoints, 0.0);
    derivCoeffs = vector<double>(npoints, 0.0);
    
    mat J = zeros<mat>(npoints,degree+1);
    for (int i = 0; i < npoints; i++) {
      for (int j = 0; j < degree+1; j++) {
	J(i,j) = pow((i-npoints/2),j);
      }
    }
    
    mat Jt = trans(J);
    mat C = solve(Jt*J,Jt);
    
    for (int i = 0; i < npoints; i++) {
      coeffs[i] = C(0,i);
      derivCoeffs[i] = C(1,i);
    }
  }
  
  ~SgolayFilter() {}

  void filter(int len, double* data, double* output) {
    if (len <= npoints) return; // Do nothing if length of data < filter size
    
    // Do filtering
    int k;
    int midpt = npoints/2;
    for (int i = 0; i < len; i++) {
      output[i] = 0.0;
      for (int j = 0; j < npoints; j++) {
	k = i+j-midpt;
	if (k < 0) {
	  k = abs(k);
	} else if (k >= len) {
	  k = len-(k-len)-1;
	}
	output[i] += data[k]*coeffs[j];
      }
    }
  }

  void derivFilter(int len, double* data, double* output) {
    if (len <= npoints) return;
    // Do filtering
    int k;
    int midpt = npoints/2;
    for (int i = 0; i < len; i++) {
      output[i] = 0.0;
      for (int j = 0; j < npoints; j++) {
	k = i+j-midpt;
	if (k < 0) {
	  k = abs(k);
	} else if (k >= len) {
	  k = len-(k-len)-1;
	}
	output[i] += data[k]*derivCoeffs[j];
      }
    }
  }
};

extern "C" {
  Filter* createSgolayFilter(int degree, int npoints) {
    return (Filter*) new SgolayFilter(degree, npoints);
  }
  
  void deleteSgolayFilter(Filter* filter) {
    delete (SgolayFilter*) filter; 
  }
  
  void sgolayFilt(Filter* filter, int len, double* data, double* output) {
    return ((SgolayFilter*) filter)->filter(len, data, output);
  }
  
  void sgolayDerivFilt(Filter* filter, int len, double* data, double* output) {
    return ((SgolayFilter*) filter)->derivFilter(len, data, output);
  }
}
