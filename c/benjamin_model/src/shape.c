// shape.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "image.h"
#include "array.h"
#include "shape.h"
#include "constant.h"

// Helper functions to compute the shape area, perimeter, and shape index
void loadData(ShapeAnalyser* ana, double** data);
void computeCM(Image* image, double* xcm, double* ycm);
void radialCoords(int npts, Point* points, double xcm, double ycm,
		  double* rad, double* theta);
void shapeGradient(int len, double* x, double* y, double* yprime);
void threePtGradient(int len, double* x, double* y, double* yprime);
int computePixelArea(int lx, int ly, double threshold, double** field);
double computeChainPerimeter(int npts, int* chain);
double computeArea(int npts, double* theta, double* r);
double computePerimeter(int npts, double* theta, double* r, double* dr);
double trapzInt(int nbins, double* x, double* y);
void insertSort(int npts, double* key, double* value);

ShapeAnalyser* createShapeAnalyser(int scale, int dataLx, int dataLy,
				   int kernelLen, double sigma,
				   int sgolayDegree, int sgolayLen, 
				   double threshold) {
  ShapeAnalyser* ana = malloc(sizeof(*ana));
  ana->scale = scale;
  ana->dataLx = dataLx;
  ana->dataLy = dataLy;
  ana->lx = dataLx*scale;
  ana->ly = dataLy*scale;
  ana->threshold = threshold;
  
  // Create the image required for the analysis
  ana->rawImage = createEmptyImage(ana->lx, ana->ly);
  ana->intermImage = createEmptyImage(ana->lx, ana->ly);
  ana->convImage = createEmptyImage(ana->lx, ana->ly);

  // Create the Gaussian blurring kernels (two 1-D kernels)
  ana->gaussX = createGaussianKernel(kernelLen, 1, sigma);
  ana->gaussY = createGaussianKernel(1, kernelLen, sigma);

  // Create Sgolay filters
  ana->sgolayDegree = sgolayDegree; 
  ana->sgolayLength = sgolayLen; // Frame length should be an odd number
  ana->sgolayRad = createSgolayFilter(ana->sgolayDegree,
				      ana->sgolayLength);
  ana->sgolayDrad = createSgolayFilter(ana->sgolayDegree,
					   ana->sgolayLength);
  return ana;
}

void deleteShapeAnalyser(ShapeAnalyser* ana) {
  deleteImage(ana->rawImage);
  deleteImage(ana->intermImage);
  deleteImage(ana->convImage);
  deleteImage(ana->gaussX);
  deleteImage(ana->gaussY);
  deleteSgolayFilter(ana->sgolayRad);
  deleteSgolayFilter(ana->sgolayDrad);
  free(ana);
}

ShapeInfo* getShapeInfo(ShapeAnalyser* ana, double** data) {
  double twopi = 2.0*PF_PI;
  
  // Rescale the phase field data to the image
  loadData(ana, data);
  
  // Blur the image using a Gaussian filter
  conv(ana->gaussX, ana->rawImage, ana->intermImage);
  conv(ana->gaussY, ana->intermImage, ana->convImage);

  // Find centre of mass
  double xcm, ycm;
  computeCM(ana->convImage, &xcm, &ycm);
  
  // Trace the boundary
  Boundary* boundary = traceBoundary(ana->threshold, ana->convImage);
  int npts = boundary->npoints; // Total number of boundary points
  int nptsp1 = npts+1;
  Point* boundpts = boundary->points;
  
  // Determine the radial coordinates
  double* theta = create1DDoubleArray(nptsp1); // +1 for endpt = starpt
  double* rad = create1DDoubleArray(nptsp1);
  radialCoords(npts, boundpts, xcm, ycm, rad, theta);
  
  // Create a periodic extension so that the smoothing process
  // can handle the boundary points as well
  int extLen = ana->sgolayLength-1; // sgolayLength should be odd
  int halfExtLen = extLen/2;
  int fullLen = extLen+npts;
  
  double* radExt = create1DDoubleArray(fullLen);
  int ext = (halfExtLen/npts+1)*npts;
  for (int i = 0; i < fullLen; i++) {
    int j = (ext+i-halfExtLen)%npts;
    radExt[i] = rad[j];
  }
  
  // Smooth radial data using Sgolay filter
  double* radExtFiltered = create1DDoubleArray(fullLen);
  sgolayFilt(ana->sgolayRad, fullLen, radExt, radExtFiltered);
  
  // Extract the non-extended part of the filtered result
  for (int i = 0; i < npts; i++) {
    rad[i] = radExtFiltered[i+halfExtLen];
  }
  
  // Add in the endpt = startpt
  theta[npts] = twopi;
  rad[npts] = rad[0];
  
  // Compute the first derivative (dr/dtheta)
  double* drad = create1DDoubleArray(nptsp1);
  shapeGradient(nptsp1, theta, rad, drad);
  drad[npts] = drad[0]; // An approximation of the gradient at the boundary
  
  // Identify numerically unstable points
  int nbadpts = 0;
  const int unstableThreshold = 100;
  double* smoothRad = create1DDoubleArray(nptsp1);
  for (int i = 0; i < npts; i++) { // Ignore endpt (which is equal to startpt)
    smoothRad[i] = rad[i];
    if (fabs(drad[i]) > unstableThreshold) {
      rad[i] = -1.0;
      nbadpts++;
    }
  } 
  
  // Do Lagrange interpolation if there are numerically unstable points
  if (nbadpts > 0) {
    int twonpts = 2*npts;
    int threenpts = 3*npts;
    const int ninter = 4;
    int nmid = ninter/2;
    double xint[ninter];
    double yint[ninter];
    double lagcoeffs[ninter];
    
    for (int i = 0; i < npts; i++) {
      // Find the unstable points
      if (rad[i] < 0.0) {
	int j = npts+i-1;
	// Find the interpolation points
	for (int k = 1; k <= nmid; k++) {
	  while (rad[j%npts] < 0.0 && j > 0.0) j--;
	  xint[nmid-k] = theta[j%npts];
	  yint[nmid-k] = rad[j%npts];
	  // Get the correct angle if the point
	  // crosses the periodic boundary
	  if (j < npts) xint[nmid-k] -= twopi;
	  j--;
	}
	j = npts+i+1;
	for (int k = nmid; k < ninter; k++) {
	  while (rad[j%npts] < 0.0 && j < threenpts-1) j++;
	  xint[k] = theta[j%npts];
	  yint[k] = rad[j%npts];
	  // Get the correct angle if the point
	  // crosses the periodic boundary
	  if (j >= twonpts) xint[k] += twopi; 
	  j++;
	}
	// Calculate the interpolation coefficients
	for (int k = 0; k < ninter; k++) {
	  lagcoeffs[k] = yint[k];
	  for (int l = 0; l < ninter; l++) {
	    if (k != l) lagcoeffs[k] /= (xint[k]-xint[l]);
	  }
	}
	// Calculate interpolated value
	smoothRad[i] = 0.0;
	for (int k = 0; k < ninter; k++) {
	  double temp = 1.0;
	  for (int l = 0; l < ninter; l++) {
	    if (k != l) temp *= (theta[i]-xint[l]);
	  }
	  smoothRad[i] += temp*lagcoeffs[k];
	}
      }
    }
    smoothRad[npts] = smoothRad[0]; // Make sure endpt = startpt
    
    // Recompute dr/dtheta
    shapeGradient(nptsp1, theta, smoothRad, drad);
    drad[npts] = drad[0];
  }

  double scale = ana->scale;
  double scale2 = scale*scale;
  ShapeInfo* info = malloc(sizeof *info);
  info->pixels = computePixelArea(ana->dataLx, ana->dataLy,
				  ana->threshold, data);
  info->area = computeArea(nptsp1, theta, smoothRad) / scale2;
  info->perimeter = computePerimeter(nptsp1, theta, smoothRad, drad) / scale;
  info->pixelArea = computePixelArea(ana->lx, ana->ly, ana->threshold,
				     ana->convImage->data) / scale2;
  info->chainPerimeter = computeChainPerimeter(npts, boundary->chain) / scale;
  
  // Clean up resources
  if (boundary != NULL) {
    deleteBoundary(boundary);
  }
  free(theta);
  free(rad);
  free(radExt);
  free(radExtFiltered);
  free(drad);
  free(smoothRad);

  return info;
}

void deleteShapeInfo(ShapeInfo* shapeInfo) {
  free(shapeInfo);
}

void loadData(ShapeAnalyser* ana, double** data) {
  for (int i = 0; i < ana->lx; i++) {
    for (int j = 0; j < ana->ly; j++) {
      ana->rawImage->data[i][j] = data[i/ana->scale][j/ana->scale];
    }
  }
}

void computeCM(Image* image, double* xcm, double* ycm) {
  *xcm = 0.0;
  *ycm = 0.0;
  double mass;
  double totalMass = 0.0;
  double** data = image->data;
  for (int i = 0; i < image->lx; i++) {
    for (int j = 0; j < image->ly; j++) {
      mass = data[i][j];
      *xcm += (i+0.5)*mass;
      *ycm += (j+0.5)*mass;
      totalMass += mass;
    }
  }
  *xcm /= totalMass;
  *ycm /= totalMass;
}

void radialCoords(int npts, Point* boundpts, double xcm, double ycm,
		  double* rad, double* theta) {
  // Get start angle
  double angle, nextAngle, startAngle, dangle;
  double pi = PF_PI;
  double twopi = 2.0*pi;
  int iangle = 0; // For counting the number of times crossed the branch cut
  double dx = boundpts[0].x+0.5-xcm;
  double dy = boundpts[0].y+0.5-ycm;
  startAngle = atan2(-dy,-dx)+pi; // Range from 0 to 2 pi
  angle = startAngle;
  theta[0] = 0.0;
  rad[0] = sqrt(dx*dx+dy*dy);

  for (int i = 1; i < npts; i++) {
    dx = boundpts[i].x+0.5-xcm;
    dy = boundpts[i].y+0.5-ycm;
    nextAngle = atan2(-dy,-dx)+pi-startAngle+twopi*iangle;
    // Check if crossed the branch cut
    dangle = nextAngle-angle;
    if (fabs(dangle) > 1.95*pi) {
      if (dangle < 0.0) {
	nextAngle += twopi;
	iangle++;
      } else {
	nextAngle -= twopi;
	iangle--;
      }
    }
    angle = nextAngle;
    theta[i] = angle;
    rad[i] = sqrt(dx*dx+dy*dy);
  }
  // Sort the angles
  insertSort(npts, theta, rad);
}

void shapeGradient(int len, double* x, double* y, double* yprime) {
  // Note that the size of yprime should be one less than that of y
  // len should be the array length of x (or y)
  // The algorithm for calculating the derivative is based on the
  // finite difference method presented in the paper "A simple
  // finite-difference grid with non-constant intervals" by Sundqvist and
  // Veronis, Tellus 22 (1970)
  int iu, id;
  int yprimeLen = len-1;
  double* dx = create1DDoubleArray(yprimeLen);
  for (int i = 0; i < yprimeLen; i++) {
    dx[i] = x[i+1]-x[i];
  }
  double dxRatio, dxRatio2;
  for (int i = 0; i < yprimeLen; i++) {
    iu = (i+1+yprimeLen)%yprimeLen;
    id = (i-1+yprimeLen)%yprimeLen;
    dxRatio = dx[i]/dx[id];
    dxRatio2 = dxRatio*dxRatio;
    yprime[i] = (y[iu]-dxRatio2*y[id]-(1-dxRatio2)*y[i])/(dx[i]*(1+dxRatio));
  }
  free(dx);
}

void threePtGradient(int len, double* x, double* y, double* yprime) {
  int i, iu, id, idd;
  int dxLen = len-1;
  double* dx = create1DDoubleArray(dxLen);
  for (i = 0; i < dxLen; i++) {
    dx[i] = x[i+1]-x[i];
  }
  double dxsum = dx[0]+dx[1];
  // Start point
  yprime[0] = -(dxsum+dx[0])/(dx[0]*dxsum)*y[0] + dxsum/(dx[0]*dx[1])*y[1] -
    dx[0]/(dx[1]*dxsum)*y[2];
  for (i = 1; i < len-1; i++) {
    iu = i+1;
    id = i-1;
    dxsum = dx[id]+dx[i];
    yprime[i] = -dx[i]/(dx[id]*dxsum)*y[id] -
      (dx[id]-dx[i])/(dx[id]*dx[i])*y[i] + dx[id]/(dx[i]*dxsum)*y[iu];
  }
  // End point
  i = len-1;
  id = len-2;
  idd = len-3;
  dxsum = dx[id]+dx[idd];
  yprime[i] = dx[id]/(dx[idd]*dxsum)*y[idd] - dxsum/(dx[idd]*dx[id])*y[id] +
    (dxsum+dx[id])/(dx[id]*dxsum)*y[i];
  free(dx);
}

int computePixelArea(int lx, int ly, double threshold, double** field) {
  int sum = 0;
  for (int i = 0; i < lx; i++)
    for (int j = 0; j < ly; j++)
      if (field[i][j] > threshold) sum++;
  return sum; 
}

double computeChainPerimeter(int npts, int* chain) {
  // Estimate the perimeter based on corner count algorithm by
  // Proffitt and Rosen (Computer Graphics and Image Processing 20(4):347 1982)
  int nnorm = 0;
  int ndiag = 0;
  int ncorner = 0;
  for (int i = 0; i < npts; i++) {
    if (chain[i] % 2 == 0) {
      nnorm++;
    } else {
      ndiag++;
    }
    if (chain[(i+1)%npts] != chain[i]) ncorner++;
  }
  return nnorm*0.980 + ndiag*1.406 - ncorner*0.091;
}

double computeArea(int npts, double* theta, double* r) {
  // Compute the area by A = 0.5 * int_0^{2pi} dtheta r^2(theta) 
  double* r2 = create1DDoubleArray(npts);
  for (int i = 0; i < npts; i++) {
    r2[i] = r[i]*r[i];
  }
  double area = fabs(trapzInt(npts, theta, r2))/2.0;
  free(r2);
  return area;
}


double computePerimeter(int npts, double* theta, double* r, double* dr) {
  // Compute the perimeter by P = int_0^{2pi} dtheta sqrt(r^2+dr^2)
  double* dl = create1DDoubleArray(npts);
  for (int i = 0; i < npts; i++) {
    dl[i] = sqrt(r[i]*r[i]+dr[i]*dr[i]);
  }
  double perimeter = fabs(trapzInt(npts, theta, dl));
  free(dl);
  return perimeter;			  
}

double trapzInt(int nbins, double* x, double* y) {
  // Numerical integration using the trapezoidal method
  double total = 0.0;
  for (int i = 0; i < nbins-1; i++) {
    total += (y[i+1]+y[i])*(x[i+1]-x[i]);
  }
  return total/2.0;
}

void insertSort(int npts, double* key, double* value) {
  // An implementation of the insertion sort algorithm for
  // sorting a key-value pair array based on the key values
  double tmpKey, tmpValue;
  int j;
  for (int i = 1; i < npts; i++) {
    tmpKey = key[i];
    tmpValue = value[i];
    j = i-1;
    while (j >= 0 && key[j] > tmpKey) {
      key[j+1] = key[j];
      value[j+1] = value[j];
      j--;
    }
    key[j+1] = tmpKey;
    value[j+1] = tmpValue;
  }
}
