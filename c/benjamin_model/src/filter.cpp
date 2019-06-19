// filter.cpp

#include "filter.h"
#include "sgolay_filter.hpp"

extern "C" {
  Filter* createSgolayFilter(int degree, int npoints) {
    return (Filter*) new SgolayFilter(degree, npoints);
  }
  void deleteSgolayFilter(Filter* filter) {
    delete (SgolayFilter*) filter;
  }
  void sgolayfilt(Filter* filter, int len, double* data, double* output) {
    return ((SgolayFilter*) filter)->filter(len, data, output);
  }
}
