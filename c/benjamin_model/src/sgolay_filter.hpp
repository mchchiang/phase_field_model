// sgolay_filter.hpp

#ifndef SGOLAY_FILTER_HPP
#define SGOLAY_FILTER_HPP

#include <vector>

class SgolayFilter {

private:
  int degree;
  int npoints;
  std::vector<double> coeffs;
  std::vector<double> derivCoeffs;

public:
  SgolayFilter(int degree, int npoints);
  ~SgolayFilter();
  void filter(int len, double* data, double* output);
  void derivFilter(int len, double* data, double* output);
};

#endif
