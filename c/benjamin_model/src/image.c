// image.c

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "image.h"
#include "array.h"

Image* createEmptyImage(int lx, int ly) {
  Image* image = malloc(sizeof *image);
  image->lx = lx;
  image->ly = ly;
  image->data = create2DDoubleArray(lx, ly);
  return image;
}

Image* createImageFromData(int lx, int ly, double** data) {
  Image* image = malloc(sizeof *image);
  image->lx = lx;
  image->ly = ly;
  image->data = create2DDoubleArray(lx, ly);
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      image->data[i][j] = data[i][j];
    }
  }
  return image;
}

void deleteImage(Image* image) {
  free(image->data);
  free(image);
}

Image* createGaussianKernel(int lx, int ly, double sigma) {
  double** data = create2DDoubleArray(lx, ly);
  int xc = lx/2;
  int yc = ly/2;
  double dx, dy;
  double sigma2 = sigma*sigma;

  // Generate the 2D Gaussian distribution
  double value; 
  double sum = 0.0;
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      dx = i-xc;
      dy = j-yc;
      value = exp(-(dx*dx+dy*dy)/(2.0*sigma2));
      data[i][j] = value;
      sum += value;
    }
  }
  
  // Normalise
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      data[i][j] /= sum;
    }
  }

  // Create the image kernel
  Image* image = malloc(sizeof *image);
  image->lx = lx;
  image->ly = ly;
  image->data = data;
  return image;
}

void conv(Image* kernel, Image* image, Image* convImage) {
  double** convImgData = convImage->data;
  double** kernelData = kernel->data;
  double** imgData = image->data;
  int xc = kernel->lx/2;
  int yc = kernel->ly/2;
  int x, y;
  
  // Do the convolution (use zero padding here for edges)
  for (int i = 0; i < image->lx; i++) {
    for (int j = 0; j < image->ly; j++) {
      convImgData[i][j] = 0.0; // Clear any previous stored data
      for (int k = 0; k < kernel->lx; k++) {
	for (int l = 0; l < kernel->ly; l++) {
	  x = i+(k-xc);
	  y = j+(l-yc);
	  if (x >= 0 && x < image->lx && y >= 0 && y < image->ly) {
	    convImgData[i][j] += imgData[x][y]*kernelData[k][l];
	  }
	}
      }
    }
  }
}

Boundary* traceBoundary(double threshold, Image* image) {
  // Retreive image data
  double** imageData = image->data;
  
  // Store the eight direction vectors
  Point dir[8] = {{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1}};

  // Array to store the boundary points
  int maxpts = image->lx*2+image->ly*2;
  Point* boundpts = malloc(sizeof *boundpts * maxpts);
  int* chain = malloc(sizeof *chain * maxpts);
  
  // Find the starting point (the first point beyond the threshold)
  Point startpt = {-1,-1};
  bool foundStart = false;
  for (int i = 0; i < image->lx && !foundStart; i++) {
    for (int j = 0; j < image->ly && !foundStart; j++) {
      if (imageData[i][j] > threshold) {
	startpt.x = i;
	startpt.y = j;
	foundStart = true;
      }
    }
  }

  if (!foundStart) {
    return NULL;
  }
  
  // Determine the first direction
  Point pt, nextpt;
  int idir = 0;
  int ptcount = 0;
  pt.x = startpt.x;
  pt.y = startpt.y;
  
  // Use the Moore-Neighbour algorithm to trace the boundary
  do {
    // Add the point to the boundary list
    boundpts[ptcount].x = pt.x;
    boundpts[ptcount].y = pt.y;
    ptcount++;
    // Extend the array if the number of boundary points exceed
    // the current length of the array
    if (ptcount >= maxpts) {
      Point* newpts;
      int* newchain;
      newpts = realloc(boundpts, sizeof *boundpts * maxpts*2);
      newchain = realloc(chain, sizeof *chain * maxpts*2);
      if (newpts == NULL || newchain == NULL) {
	if (newpts != NULL) free(newpts);
	if (newchain != NULL) free(newchain);
	break;
      }
      boundpts = newpts;
      chain = newchain;
      maxpts *= 2;
    }
    // Search for the next boundary point in an anti-clockwise direction
    idir -= 4; // Direction to backtrack to previous point
    if (idir < 0) idir += 8;
    else if (idir >= 8) idir -= 8;
    for (int i = 0; i < 8; i++) {
      idir = (idir+1)%8;
      nextpt.x = pt.x+dir[idir].x;
      nextpt.y = pt.y+dir[idir].y;
      if (nextpt.x >= 0 && nextpt.x < image->lx &&
	  nextpt.y >= 0 && nextpt.y < image->ly &&
	  imageData[nextpt.x][nextpt.y] > threshold) {
	// Found the next boundary point
	pt.x = nextpt.x;
	pt.y = nextpt.y;
	chain[ptcount-1] = idir; 
	break;
      }
    }
  }  while (pt.x != startpt.x || pt.y != startpt.y);

  // Resize the array of boundary points to the actual size of the array
  Point* newpts;
  int* newchain;
  newpts = realloc(boundpts, sizeof *boundpts * ptcount);
  newchain = realloc(chain, sizeof *chain * ptcount);
  if (newpts == NULL || newchain == NULL) {
    if (newpts != NULL) free(newpts);
    if (newchain != NULL) free(newchain);
    Boundary* boundary = malloc(sizeof *boundary);
    boundary->points = boundpts;
    boundary->chain = chain;
    boundary->npoints = ptcount;
    return boundary;
  }
  boundpts = newpts;
  chain = newchain;
  
  Boundary* boundary = malloc(sizeof *boundary);
  boundary->points = boundpts;
  boundary->chain = chain;
  boundary->npoints = ptcount;
  return boundary;
}

void deleteBoundary(Boundary* boundary) {
  free(boundary->points);
  free(boundary->chain);
  free(boundary);
}
