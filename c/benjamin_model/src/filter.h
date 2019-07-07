// filter.h

#ifndef FILTER_H
#define FILTER_H

#ifdef __cplusplus
extern "C" {
#endif
  
#include <stdio.h>
#include <stdlib.h>

  typedef struct Filter Filter;
  Filter* createSgolayFilter(int degree, int npoints);
  void deleteSgolayFilter(Filter* filter);
  void sgolayFilt(Filter* filter, int len, double* data, double* output);
  void sgolayDerivFilt(Filter* filter, int len, double* data, double* output);
  
#ifdef __cplusplus
}
#endif
#endif
