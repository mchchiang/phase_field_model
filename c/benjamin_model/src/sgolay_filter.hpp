// sgolay_filter.hpp

#ifndef SGOLAY_FILTER_HPP
#define SGOLAY_FILTER_HPP

#include <vector>

struct SgolayFilter {
  int degree;
  int npoints;
  std::vector<double> coeffs;

  SgolayFilter(int degree, int npoints);
  ~SgolayFilter();
  void filter(int len, double* data, double* output);
};

#endif
