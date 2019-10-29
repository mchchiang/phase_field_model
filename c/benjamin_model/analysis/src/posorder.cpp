// posorder.cpp
// A program to compute the positional order parameter

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
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
template <typename T> int sgn(T val);

int main(int argc, char* argv[]) {
  
  if (argc != 9) {
    cout << "usage: posorder npoints lx ly startTime endTime timeInc "
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

  // First reciprocal lattice vector
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }

  const double R {8.0};
  const double twopi {2.0*M_PI};
  const double qx {twopi/R};
  const double qy {twopi/R*sqrt(3.0)};
  
  double x, y, dx, dy, sum;
  double nsq {static_cast<double>(npoints*npoints)};
  
  for (ibin = 0; ibin < nbins; ibin++) {
    sum = 0.0;
    for (int i {}; i < npoints; i++) {
      x = (*pos)[ibin][i][0];
      y = (*pos)[ibin][i][1];	
      for (int j {}; j < i; j++) {
	dx = diff(lx, x, (*pos)[ibin][j][0]);
	dy = diff(ly, y, (*pos)[ibin][j][1]);
	sum += 2.0*cos(dot(qx, qy, dx, dy));
      }
    }
    sum = (sum + npoints)/nsq;
    time = ibin*timeInc+startTime;
    writer << time << " " << sum << endl;
  }
  
  // Clean up resources
  delete pos;
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
