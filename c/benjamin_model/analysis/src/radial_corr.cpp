// radial_corr.cpp
// A program to compute the radial correlation function

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
    cout << "usage: radial_corr npoints lx ly startTime endTime timeInc "
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

  double x, y, dx, dy, dr;
  vector<double> nvec (npoints, 0.0);
  vector<vector<double> >* drMat {new vector<vector<double> >(npoints, nvec)};
  
  // Compute relative position matrix
  for (int i {}; i < npoints; i++) {
    x = (*pos)[ibin][i][0];
    y = (*pos)[ibin][i][1];	
    for (int j {}; j < i; j++) {
      dx = diff(lx, x, (*pos)[ibin][j][0]);
      dy = diff(ly, y, (*pos)[ibin][j][1]);
      dr = dist(dx, dy);
      (*drMat)[i][j] = dr;
      (*drMat)[j][i] = dr;
    }
  }
  const double twopi {2.0*M_PI};
  const double incr {0.1};
  const double rmax {lx<ly?lx/2.0:ly/2.0};
  int nr {static_cast<int>(floor(rmax/incr))};
  double r1, r2, sum {};
  for (int ir {}; ir < nr; ir++) {
    sum = 0.0;
    r1 = ir*incr;
    r2 = r1+incr;
    for (int i {}; i < npoints; i++) {
      for (int j {}; j < npoints; j++) {
	if (i == j) continue; 
	if ((*drMat)[i][j] >= r1 && (*drMat)[i][j] < r2) {
	  sum += 1.0;
	}
      }
    }
    sum /= (static_cast<int>(npoints)*twopi*r1*incr);
    writer << r1 << " " << r2 << " " << 0.5*(r1+r2) << " " << sum << endl;
  }

  // Clean up resources
  delete pos;
  delete drMat;
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
