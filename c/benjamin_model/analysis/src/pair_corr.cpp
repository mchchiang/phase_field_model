// pair_corr.cpp
// A program to compute the pair correlation function

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <omp.h>
#include "position.hpp"

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
double gauss(double sigma, double x, double y);
template <typename T> int sgn(T val);

int main(int argc, char* argv[]) {
  
  if (argc != 9) {
    cout << "usage: structure npoints lx ly startTime endTime timeInc "
	 << "posFile outFile" << endl;
    return 1;
  }

  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
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
  vector<double> vec (2, 0.0);
  vector<vector<double> > vec2 (npoints, vec);
  vector<vector<vector<double> > >* pos 
  {new vector<vector<vector<double > > >(nbins, vec2)};
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
	(*pos)[ibin][i][0] = posReader.getPosition(i, 0);
	(*pos)[ibin][i][1] = posReader.getPosition(i, 1);
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

  ibin = 0;

  double x, y;
  vector<vector<vector<double> > >* diffMat 
  {new vector<vector<vector<double> > >(npoints, vec2)};
  
  // Compute relative position matrix
  for (int i {}; i < npoints; i++) {
    x = (*pos)[ibin][i][0];
    y = (*pos)[ibin][i][1];	
    for (int j {}; j < i; j++) {
      (*diffMat)[i][j][0] = diff(lx, x, (*pos)[ibin][j][0]);
      (*diffMat)[i][j][1] = diff(ly, y, (*pos)[ibin][j][1]);
      (*diffMat)[j][i][0] = -(*diffMat)[i][j][0];
      (*diffMat)[j][i][1] = -(*diffMat)[i][j][1];
    }
  }
  
  const double sigma {0.01};
  const double d {0.1};
  const int nx {static_cast<int>(ceil(lx/d))};
  const int ny {static_cast<int>(ceil(ly/d))};
  double dx, dy, dx0, dy0, sum;
  for (int i {}; i < nx; i++) {
    dx = i*d-lx/2.0;
    for (int j {}; j < ny; j++) {
      dy = j*d-ly/2.0;
      sum = 0.0;
      cout << i << " " << j << endl;
      for (int k {}; k < npoints; k++) {
	for (int l {}; l < k; l++) {
	  dx0 = (*diffMat)[k][l][0];
	  dy0 = (*diffMat)[k][l][1];
	  sum += gauss(sigma, (dx-dx0), (dy-dy0));
	}
      }
      writer << dx << " " << dy << " " << sum << endl;
    }
  }

  // Clean up resources
  delete pos;
  delete diffMat;
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

double gauss(double sigma, double x, double y) {
  return exp(-(x*x+y*y)/(sigma*sigma));
}

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
