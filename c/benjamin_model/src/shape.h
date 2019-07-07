// shape.h
// Codes for analysing the geometry/shape of a phase field

#ifndef SHAPE_H
#define SHAPE_H

#include "image.h"
#include "filter.h"

typedef struct ShapeAnalyser {
  int scale;
  int dataLx;
  int dataLy;
  int lx;
  int ly;
  double threshold;

  // For image analysis
  Image* rawImage;
  Image* intermImage;
  Image* convImage;

  // Gaussian blurring kernels
  Image* gaussX;
  Image* gaussY;

  // For signal filtering
  int sgolayDegree;
  int sgolayLength;
  Filter* sgolayRad;
  Filter* sgolayDrad;
} ShapeAnalyser;

typedef struct ShapeInfo {
  int pixels; 
  double area;
  double pixelArea;
  double perimeter;
  double chainPerimeter;
} ShapeInfo;

ShapeAnalyser* createShapeAnalyser(int scale, int dataLx, int dataLy,
				   int kernelLength, double sigma,
				   int sgolayDegree, int sgolayLength, 
				   double threshold);
void deleteShapeAnalyser(ShapeAnalyser* ana);
ShapeInfo* getShapeInfo(ShapeAnalyser* ana, double** data);
void deleteShapeInfo(ShapeInfo* shapeInfo);

#endif
