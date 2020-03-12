// array.cpp

#include "arralloc.h"
#include "array.hpp"

double* create1DDoubleArray(int len) {
  double* arr = (double*) arralloc(sizeof(double), 1, len);
  for (int i = 0; i < len; i++) {
    arr[i] = 0.0;
  }
  return arr;
}

double** create2DDoubleArray(int len1, int len2) {
  double** arr = (double**) arralloc(sizeof(double), 2, len1, len2);
  for (int i = 0; i < len1; i++) {
    for (int j = 0; j < len2; j++) {
      arr[i][j] = 0.0;
    }
  }
  return arr;
}

double*** create3DDoubleArray(int len1, int len2, int len3) {
  double*** arr = (double***) arralloc(sizeof(double), 3, len1, len2, len3);
  for (int i = 0; i < len1; i++) {
    for (int j = 0; j < len2; j++) {
      for (int k = 0; k < len3; k++) {
	arr[i][j][k] = 0.0;
      }
    }
  }
  return arr;
}

int* create1DIntArray(int len) {
  int* arr = (int*) arralloc(sizeof(int), 1, len);
  for (int i = 0; i < len; i++) {
    arr[i] = 0;
  }
  return arr;
}

int** create2DIntArray(int len1, int len2) {
  int** arr = (int**) arralloc(sizeof(int), 2, len1, len2);
  for (int i = 0; i < len1; i++) {
    for (int j = 0; j < len2; j++) {
      arr[i][j] = 0;
    }
  }
  return arr;
}

int*** create3DIntArray(int len1, int len2, int len3) {
  int*** arr = (int***) arralloc(sizeof(int), 3, len1, len2, len3);
  for (int i = 0; i < len1; i++) {
    for (int j = 0; j < len2; j++) {
      for (int k = 0; k < len3; k++) {
	arr[i][j][k] = 0;
      }
    }
  }
  return arr;
}
