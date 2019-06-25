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
double computeArea(int npts, double* theta, double* r);
double computePerimeter(int npts, double* theta, double* r, double* dr);
double trapzInt(int nbins, double* x, double* y);
void insertSort(int npts, double* key, double* value);

ShapeAnalyser* createShapeAnalyser(int scale, int dataLx, int dataLy,
				   int kernelLen, double sigma,
				   int sgolayDegree, int sgolayLen, 
				   double threshold) {
  ShapeAnalyser* ana = (ShapeAnalyser*) malloc(sizeof(*ana));
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
  /*double sigma = 4.0;
  int sigmaInt = (int) (ceil(sigma));
  int kernelL = 6*sigmaInt-sigmaInt%2+1;*/
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

void getShapeInfo(ShapeAnalyser* ana, double** data, double* perimeter,
		  double* area, double* asphericity) {
  // Rescale the phase field data to the image
  loadData(ana, data);

  /*FILE* f = fopen("raw_c.dat", "w");
  for (int i = 0; i < ana->lx; i++) {
    for (int j = 0; j < ana->ly; j++) {
      fprintf(f, "%d %d %.5f\n", i, j, ana->rawImage->data[i][j]);
    }
  }
  fclose(f);*/
  
  // Blur the image using a Gaussian filter
  conv(ana->gaussX, ana->rawImage, ana->intermImage);
  conv(ana->gaussY, ana->intermImage, ana->convImage);

  // Output the convolved image
  /*f = fopen("conv_c.dat", "w");
  for (int i = 0; i < ana->lx; i++) {
    for (int j = 0; j < ana->ly; j++) {
      fprintf(f, "%d %d %.5f\n", i, j, ana->convImage->data[i][j]);
    }
  }
  fclose(f);*/
  
  // Find centre of mass
  double xcm, ycm;
  computeCM(ana->convImage, &xcm, &ycm);
    
  // Trace the boundary
  int npts; // Total number of boundary points
  Point* boundpts = traceBoundary(ana->threshold, ana->convImage, &npts);
  
  // Output boundary points
  /*f = fopen("pts_c.dat", "w");
  for (int i = 0; i < npts; i++) {
    fprintf(f, "%d %d\n", boundpts[i].x, boundpts[i].y);
  }
  fclose(f);*/
  
  // Determine the radial coordinates
  double* theta = create1DDoubleArray(npts);
  double* rad = create1DDoubleArray(npts);
  radialCoords(npts, boundpts, xcm, ycm, rad, theta);
  
  // Create a periodic extension so that the smoothing process
  // can handle the end points well
  int extLen = ana->sgolayLength-1;
  //int filtLen = ((npts*2-1)/2)*2+1; // Make sure it is an odd number
  //int extLen = filtLen-1;
  int halfExtLen = extLen/2;
  int fullLen = extLen+npts;
  double* radExt = create1DDoubleArray(fullLen);
  for (int i = 0; i < fullLen; i++) {
    int j = (npts+i-halfExtLen)%npts;
    radExt[i] = rad[j];
  }
  
  // Smooth radial data using Sgolay filter
  double* radExtFiltered = create1DDoubleArray(fullLen);
  //Filter* sgofilt = createSgolayFilter(ana->sgolayDegree, filtLen);  
  sgolayfilt(ana->sgolayRad, fullLen, radExt, radExtFiltered);
  //sgolayfilt(sgofilt, fullLen, radExt, radExtFiltered);

  // Output the raw and smoothed radial data
  /*f = fopen("radial_c.dat", "w");
  for (int i = 0; i < npts; i++) {
    fprintf(f, "%.5f %.5f %.5f\n", theta[i], rad[i],
	    radExtFiltered[i+halfExtLen]);
  }
  fclose(f);*/
  
  // Extract the non-extended part of the filtered result
  for (int i = 0; i < npts; i++) {
    rad[i] = radExtFiltered[i+halfExtLen];
  }

  // Compute the first derivative (dr/dtheta)
  int nptsm1 = npts-1;
  double* drad = create1DDoubleArray(nptsm1);
  shapeGradient(npts, theta, rad, drad);
  
  // Smooth the gradient data
  for (int i = 0; i < nptsm1; i++) {
    int iu = (i+1+nptsm1)%(nptsm1);
    int id = (i-1+nptsm1)%(nptsm1);
    if (fabs(drad[i]) > 100) {
      drad[i] = (drad[iu]+drad[id])/2.0;
    }
  }
  
  double* dradExt = create1DDoubleArray(fullLen-1);
  double* dradExtFiltered = create1DDoubleArray(fullLen-1);
  for (int i = 0; i < fullLen-1; i++) {
    int j = (npts+i-halfExtLen-1)%(nptsm1);
    dradExt[i] = drad[j];
  }
  
  sgolayfilt(ana->sgolayDrad, fullLen-1, dradExt, dradExtFiltered);
  //sgolayfilt(sgofilt, fullLen-1, dradExt, dradExtFiltered);
  //deleteSgolayFilter(sgofilt);
  // Output the smoothed radial data
  /*f = fopen("grad_radial_c.dat", "w");
  for (int i = 0; i < nptsm1; i++) {
    fprintf(f, "%.5f %.5f %.5f\n", theta[i], drad[i],
	    dradExtFiltered[i+halfExtLen]);
  }
  fclose(f);*/

  for (int i = 0; i < nptsm1; i++) {
    drad[i] = dradExtFiltered[i+halfExtLen];
  }

  // Compute area
  *area = computeArea(nptsm1, theta, rad)/(ana->scale*ana->scale);
  
  // Compute perimeter
  *perimeter = computePerimeter(nptsm1, theta, rad, drad)/ana->scale;
  
  // Calculate the perimeter to area ratio
  *asphericity = (*perimeter)*(*perimeter)/((*area)*4.0*PF_PI);
  
  // Clean up resources
  if (boundpts != NULL) {
    free(boundpts);
  }
  free(theta);
  free(rad);
  free(radExt);
  free(radExtFiltered);
  free(dradExtFiltered);
  free(drad);
  free(dradExt);
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
    dxRatio = dx[i]-dx[id];
    dxRatio2 = dxRatio*dxRatio;
    yprime[i] = (y[iu]-dxRatio2*y[id]-(1-dxRatio2)*y[i])/(dx[i]*(1+dxRatio));
  }
  free(dx);
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
    total += (y[i+1]+y[i])*(x[i+1]-x[i])/2.0;
  }
  return total;
}

// An implementation of the insertion sort algorithm for
// sorting a key-value pair array based on the key values
void insertSort(int npts, double* key, double* value) {
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
