// sgolay_filter.cpp

#include <cmath>
#include <vector>
#include <armadillo>
#include "sgolay_filter.hpp"

using std::vector;
using namespace arma;

SgolayFilter::SgolayFilter(int z, int npts) {
  degree = z;
  npoints = npts;

  // Find the coefficients
  coeffs = vector<double>(npoints, 0.0);

  mat J = zeros<mat>(npoints,degree+1);
  for (int i = 0; i < npoints; i++) {
    for (int j = 0; j < degree+1; j++) {
      J(i,j) = pow((i-npoints/2),j);
    }
  }

  mat Jt = trans(J);
  mat JtJinv = inv(Jt*J);
  mat C = JtJinv*Jt;

  for (int i = 0; i < npoints; i++) {
    coeffs[i] = C(0,i);
  }
}

SgolayFilter::~SgolayFilter() {}

void SgolayFilter::filter(int len, double* data, double* output) {
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
