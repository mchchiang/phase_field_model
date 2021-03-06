// structure.cpp
// A program to compute the structure factor

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <random>
#include <omp.h>
#include "position.hpp"
#include "array.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::string;
using std::vector;

double diff(double lx, double x1, double x2);
double dist(double x, double y);
double dot(double x1, double y1, double x2, double y2);
template <typename T> int sgn(T val);

int main(int argc, char* argv[]) {
  
  if (argc != 10) {
    cout << "usage: structure npoints nqvec lx ly startTime endTime timeInc "
	 << "posFile outFile" << endl;
    return 1;
  }

  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
  int nqvec {stoi(string(argv[++argi]), nullptr, 10)};
  int lx {stoi(string(argv[++argi]), nullptr, 10)};
  int ly {stoi(string(argv[++argi]), nullptr, 10)};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  long timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  string posFile {argv[++argi]};
  string outFile {argv[++argi]};
  
  PositionReader posReader;
  posReader.open(posFile, npoints, lx, ly, timeInc);
  if (!posReader.isOpen()) {
    cout << "Problem with opening position file!" << endl;
    return 1;
  }

  // Read the position data
  int nbins {static_cast<int>((endTime-startTime)/timeInc)+1};
  double*** pos {create3DDoubleArray(nbins, npoints, 2)};
  long time;
  int ibin;
  while (posReader.nextFrame()) {
    time = posReader.getTime();
    if (time < startTime) {
      continue;
    } else if (time >= startTime && time <= endTime) {
      ibin = static_cast<int>((time-startTime)/timeInc);
      for (int i {}; i < npoints; i++) {
	// Use wrapped position here as we only care about relative distances
	// and orientations
	pos[ibin][i][0] = posReader.getPosition(i, 0);
	pos[ibin][i][1] = posReader.getPosition(i, 1);
      }
    } else {
      break;
    }
  }
  posReader.close();
  
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening average output file!" << endl;
    return 1;
  }

  // Computing the structure factor
  double twopi {2.0*M_PI};
  double R {8.0};
  double qxmin {-twopi/R};
  double qxmax {twopi/R};
  double qymin {-twopi/R};
  double qymax {twopi/R};
  double qlx {qxmax-qxmin};
  double qly {qymax-qymin};
  double qxinc {qlx/static_cast<double>(nqvec)};
  double qyinc {qly/static_cast<double>(nqvec)};
  double qx, qy, sq;

  double*** diffMat {create3DDoubleArray(npoints, npoints, 2)};
  double** structFact {create2DDoubleArray(nqvec, nqvec)};
  double x, y, prod;
  
  for (ibin = 0; ibin < nbins; ibin++) {
    cout << "Doing t = " << (startTime+timeInc*ibin) << endl;
    // Compute the separation between points
    for (int i = 0; i < npoints; i++) {
      x = pos[ibin][i][0];
      y = pos[ibin][i][1];	
      for (int j = 0; j < i; j++) {
	diffMat[i][j][0] = diff(lx, x, pos[ibin][j][0]);
	diffMat[i][j][1] = diff(ly, y, pos[ibin][j][1]);
      }
    }
      
#pragma omp parallel default(none),\
  shared(diffMat, qxmin, qymin, qxinc, qyinc, nqvec, npoints, structFact),\
  private(qx, qy, sq, prod)
    {
#pragma omp for schedule(static) collapse(2)
      for (int i = 0; i < nqvec; i++) {
	for (int j = 0; j < nqvec; j++) {
	  qx = qxmin+qxinc*(i+0.5);
	  qy = qymin+qyinc*(j+0.5);
	  sq = 0.0;
	  for (int k = 0; k < npoints; k++) {
	    for (int l = 0; l < k; l++) {
	      prod = dot(qx, qy, diffMat[k][l][0], diffMat[k][l][1]);
	      sq += cos(prod);
	    }
	  }
	  sq = sq*2.0/npoints+1.0;
	  structFact[i][j] += sq;
	}
      }
    }
  }
  
  // Normalise and output results
  writer << std::setprecision(10) << std::fixed;
  for (int i = 0; i < nqvec; i++) {
    for (int j = 0; j < nqvec; j++) {
      structFact[i][j] /= static_cast<double>(nbins);
      qx = qxmin+qxinc*(i+0.5);
      qy = qymin+qyinc*(j+0.5);
      writer << qx << " " << qy << " " << structFact[i][j] << "\n";
    }
    writer << "\n"; 
  }
  writer.close();

  // Clean up resources
  free(pos);
  free(diffMat);
  free(structFact);
}

double diff(double lx, double x1, double x2) {
  double dx1 {x1-x2};
  double dx1mag {fabs(dx1)};
  double dx2mag {lx-dx1mag};
  return dx2mag > dx1mag ? dx1 : -dx2mag*sgn(dx1);
}

double dist(double x, double y) {
  return sqrt(x*x+y*y);
}

double dot(double x1, double y1, double x2, double y2) {
  return x1*x2+y1*y2;
}

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
