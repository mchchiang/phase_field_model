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

Image* createEmptyImage(int lx, int ly);
Image* createImageFromData(int lx, int ly, double** data);
void deleteImage(Image* image);
Image* createGaussianKernel(int lx, int ly, double sigma);
void conv(Image* image1, Image* image2, Image* convImage);
Point* traceBoundary(double threshold, Image* image, int* npts);

#endif
