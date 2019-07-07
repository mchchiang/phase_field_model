// image.h
// Functions for image analysis

#ifndef IMAGE_H
#define IMAGE_H

typedef struct Image {
  int lx;
  int ly;
  double** data;
} Image;

typedef struct Point {
  int x;
  int y;
} Point;

typedef struct Boundary {
  int npoints; // Number of points
  Point* points; // The boundary points
  int* chain; // The chain code
} Boundary; 

Image* createEmptyImage(int lx, int ly);
Image* createImageFromData(int lx, int ly, double** data);
Image* createGaussianKernel(int lx, int ly, double sigma);
void deleteImage(Image* image);

Boundary* traceBoundary(double threshold, Image* image);
void deleteBoundary(Boundary* boundary);

void conv(Image* image1, Image* image2, Image* convImage);

#endif
