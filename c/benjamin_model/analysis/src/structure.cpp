// structure.cpp
// A program to compute the structure factor

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

  // Computing the structure factor
  const int nqvec {1000};
  double twopi {2.0*M_PI};
  double R {8.0};
  double qxmin {-twopi/R};
  double qxmax {twopi/R};
  double qymin {-twopi/R};
  double qymax {twopi/R};
  double qlx {qxmax-qxmin};
  double qly {qymax-qymin};
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> rand(0.0, 1.0);
  vector<vector<double> > diffVec (npoints, vec);
  vector<vector<vector<double> > >* diffMat
  {new vector<vector<vector<double> > >(npoints, diffVec)};
  vector<vector<double> > qvec (nqvec, vec);
  vector<double> structFact (nqvec, 0.0);
  double r, x, y, prod;
  int i, j, k;
  for (ibin = 0; ibin < nbins; ibin++) {
    // cout << "Doing bin " << ibin << endl;
    // Generate the random q vectors
    for (i = 0; i < nqvec; i++) {
      r = rand(mt);
      qvec[i][0] = qxmin+qlx*r;
      r = rand(mt);
      qvec[i][1] = qymin+qly*r;
    }

    for (i = 0; i < npoints; i++) {
      x = (*pos)[ibin][i][0];
      y = (*pos)[ibin][i][1];	
      for (j = 0; j < i; j++) {
	(*diffMat)[i][j][0] = diff(lx, x, (*pos)[ibin][j][0]);
	(*diffMat)[i][j][1] = diff(ly, y, (*pos)[ibin][j][1]);
      }
    }
      
#pragma omp parallel default(none),			\
  shared(diffMat, qvec, structFact, npoints),		\
  private(i, j, k, prod)
    {
#pragma omp for schedule(static)
      for (i = 0; i < nqvec; i++) {
	structFact[i] = 0.0;
	for (j = 0; j < npoints; j++) {
	  for (k = 0; k < j; k++) {
	    prod = dot(qvec[i][0], qvec[i][1], 
		       (*diffMat)[j][k][0], (*diffMat)[j][k][1]);
	    structFact[i] += cos(prod);
	  }
	}
	structFact[i] = structFact[i]*2.0/npoints+1.0;
      }
    }
    
    // Output the points
    for (i = 0; i < nqvec; i++) {
      writer << qvec[i][0] << " " << qvec[i][1] << " " 
	     << structFact[i] << endl;
    }
  }
  writer.close();

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

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
