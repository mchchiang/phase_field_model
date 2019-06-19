// test_image.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "array.h"
#include "image.h"
#include "filter.h"

double trapzInt(int nbins, double* x, double* y);

int main(int argc, char* argv[]) {
  // Make a circle
  int lx = 164;
  int ly = 164;
  int xc = lx/2;
  int yc = ly/2;
  double dx, dy;
  double r2 = 32*32;
  Image* image = createEmptyImage(lx, ly);
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      dx = i-xc;
      dy = j-yc;
      if (dx*dx+dy*dy < r2) {
	image->data[i][j] = 2.0;
      }
    }
  }

  FILE* f = fopen("circle_original.dat", "w");
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      fprintf(f, "%d %d %.5f\n", i, j, image->data[i][j]);
    }
    fprintf(f,"\n");
  }
  fclose(f);
  
  // Create the Gaussian kernel for blurring
  int kernelLx = 23;
  int kernelLy = 23;
  double sigma = 4.0;
  Image* gaussX = createGaussianKernel(kernelLx, 1, sigma);
  Image* gaussY = createGaussianKernel(1, kernelLy, sigma);
  
  // Do image convolution
  Image* convImg = createEmptyImage(lx, ly);
  Image* intImg = createEmptyImage(lx, ly);
  
  conv(gaussX, image, intImg);
  conv(gaussY, intImg, convImg);
  
  f = fopen("circle_gauss.dat", "w");
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      fprintf(f, "%d %d %.5f\n", i, j, convImg->data[i][j]);
    }
    fprintf(f,"\n");
  }
  fclose(f);
  
  // Trace the boundary
  int npts = 0;
  Point* boundpts = traceBoundary(1.0, convImg, &npts);
  printf("%d\n", npts);
  f = fopen("bound_pts.dat", "w");
  
  for (int i = 0; i < npts; i++) {
    fprintf(f, "%d %d\n", boundpts[i].x, boundpts[i].y);
  }
  fclose(f);
  printf("Done\n");
  
  // Find centre of mass
  double xcm = 0.0;
  double ycm = 0.0;
  double mass;
  double totalMass = 0.0;
  double** convImgData = convImg->data;
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      mass = convImgData[i][j];
      if (mass > 1.0) {
	xcm += (i+0.5)*mass;
	ycm += (j+0.5)*mass;
	totalMass += mass;
      }
    }
  }
  xcm /= totalMass;
  ycm /= totalMass;
  
  // Compute radial coordinates
  double* theta = create1DDoubleArray(npts);
  double* rad = create1DDoubleArray(npts);

  // Get start angle
  double angle, nextAngle, startAngle, dangle;
  double pi = M_PI;
  double twopi = 2.0*M_PI;
  int iangle = 0; // For counting the number of times crossed the branch cut
  dx = boundpts[0].x+0.5-xcm;
  dy = boundpts[0].y+0.5-ycm;
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
    if (fabs(dangle) > 1.95*M_PI) {
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

  for (int i = 0; i < npts; i++) {
    printf("%d %.5f %.5f\n", i, theta[i], rad[i]);
  }
  
  // Smooth the radial data
  Filter* sgolay = createSgolayFilter(5,11);
  double* radf = create1DDoubleArray(npts);
  /*double* radext = create1DDoubleArray(npts*2);
  for (int i = 0; i < npts*2; i++) {
    radext[i] = rad[(i+npts/2)%npts];
    }*/
  sgolayfilt(sgolay, npts, rad, radf);
  
  for (int i = 0; i < npts; i++) {
    printf("%d %.5f %.5f\n", i, theta[i], radf[i]);
  }
  
  double radAvg = 0.0;
  for (int i = 0; i < npts; i++) {
    radAvg += radf[i];
  }
  radAvg /= (double) npts;
  printf("radial average: %.10f\n",radAvg);
  
  // Estimate dr/dtheta
  double* dtheta = create1DDoubleArray(npts);
  double* drdtheta = create1DDoubleArray(npts);
  double* drdthetaf = create1DDoubleArray(npts);
  for (int i = 0; i < npts; i++) {
    int j = (i+1)%npts;
    dtheta[i] = theta[j]-theta[i];
  }
  
  // Use the method proposed in (Ref.) to estimate the derivative
  double thetaRatio, thetaRatio2;
  for (int i = 0; i < npts; i++) {
    int iu = (i+1+npts)%npts;
    int id = (i-1+npts)%npts;
    thetaRatio = dtheta[i]-dtheta[id];
    thetaRatio2 = thetaRatio*thetaRatio;
    drdtheta[i] = (radf[iu]-thetaRatio2*radf[id]-(1.0-thetaRatio2)*radf[i]) /
      (dtheta[i]*(1+thetaRatio));
  }
  sgolayfilt(sgolay, npts, drdtheta, drdthetaf);

  // Compute area
  double* rad2 = create1DDoubleArray(npts);
  for (int i = 0; i < npts; i++) {
	rad2[i] = radf[i]*radf[i];
  }
  double area = fabs(trapzInt(npts, theta, rad2)) / 2.0;
  printf("area: %.10f\n", area);

  // Compute perimeter
  double* dl = create1DDoubleArray(npts);
  for (int i = 0; i < npts; i++) {
	dl[i] = sqrt(radf[i]*radf[i]+drdtheta[i]*drdtheta[i]);
  }
  double perimeter = fabs(trapzInt(npts, theta, dl));
  printf("perimeter: %.10f\n", perimeter);
  printf("ratio: %.10f\n", perimeter/sqrt(area));
  
  // Clean up resources
  free(dtheta);
  free(drdtheta);
  free(drdthetaf);
  free(rad);
  free(theta);
  free(rad2);
  free(dl);
  free(radf);
  if (boundpts != NULL) {
    free(boundpts);
  }
  deleteImage(gaussX);
  deleteImage(gaussY);
  deleteImage(image);
  deleteImage(intImg);
  deleteImage(convImg);
  deleteSgolayFilter(sgolay);
}

double trapzInt(int nbins, double* x, double* y) {
  double total = 0.0;
  for (int i = 0; i < nbins-1; i++) {
    total += (y[i+1]+y[i])*(x[i+1]-x[i])/2.0;
  }
  return total;
}
