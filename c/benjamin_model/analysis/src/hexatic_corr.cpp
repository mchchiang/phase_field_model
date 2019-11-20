// hexatic_corr.cpp
// A program to compute the correlation of the hexatic order parameter

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <set>
#include <utility>
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
using std::set;
using std::pair;

double diff(double lx, double x1, double x2);
double dist(double x, double y);
double dot(double x1, double y1, double x2, double y2);
template <typename T> int sgn(T val);

int main(int argc, char* argv[]) {
  
  if (argc != 13) {
    cout << "usage: hexatic_corr npoints lx ly min max binSize startTime "
	 << "endTime timeInc posFile neighFile outFile" << endl;
    return 1;
  }

  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
  int lx {stoi(string(argv[++argi]), nullptr, 10)};
  int ly {stoi(string(argv[++argi]), nullptr, 10)};
  double min {stod(string(argv[++argi]), nullptr)};
  double max {stod(string(argv[++argi]), nullptr)};
  double binSize {stod(string(argv[++argi]), nullptr)};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  long timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  string posFile {argv[++argi]};
  string neighFile {argv[++argi]};
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
	(*pos)[ibin][i][0] = posReader.getPosition(i, 0);
	(*pos)[ibin][i][1] = posReader.getPosition(i, 1);
      }
    } else {
      break;
    }
  }
  posReader.close();

  // Read the neighbour data
  vector<set<pair<int,int> > >* bonds 
  {new vector<set<pair<int,int> > >(nbins)};
  ifstream neighReader;
  neighReader.open(neighFile);
  if (!neighReader) {
    cout << "Problem with opening neighbour file!" << endl;
    return 1;
  }
  string line, timestr;
  time = 0L;
  istringstream iss;
  int neighIndex;
  while (getline(neighReader, line)) {
    // Read time
    if (getline(neighReader, line)) {
      iss.clear();
      iss.str(line);
      iss >> timestr >> time;      
      if (time < startTime) {
	// Skip over unused data
	for (int i {}; i < npoints; i++) {
	  if (!getline(neighReader, line)) break;
	}
      } else if (time >= startTime && time <= endTime) {
	ibin = static_cast<int>((time-startTime)/timeInc);
	for (int i {}; i < npoints; i++) {
	  getline(neighReader, line);
	  iss.clear();
	  iss.str(line);
	  while (iss) {
	    iss >> neighIndex;
	    if (!iss) break;
	    if (neighIndex > i) {
	      (*bonds)[ibin].insert(std::make_pair(i, neighIndex));
	    } else {
	      (*bonds)[ibin].insert(std::make_pair(neighIndex, i));
	    }
	  }
	}
      } else {
	break;
      }
    } else {
      break;
    }
  }
  neighReader.close();
  
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  writer << std::setprecision(5) << std::fixed;
  
  int idbin;
  int ndbins {static_cast<int>(ceil((max-min)/binSize))};
  vector<vector<double> > hexDistrb(ndbins, vec);
  vector<int> distrbCount(ndbins, 0);
  int ind1, ind2;
  int nbonds;
  double xi, xj, yi, yj, dx, dy, dr, theta;
  double hexxi, hexxj, hexyi, hexyj, corrx, corry;
  const double pi {M_PI};
  
  // Compute hexatic correlation
  for (ibin = 0; ibin < nbins; ibin++) {
    int n {};
    // Compute bond orientational order parameter
    // Check the actual number of bonds
    nbonds = (*bonds)[ibin].size();
    vector<vector<double> >* posBond 
    {new vector<vector<double> >(nbonds,vec)};
    vector<vector<double> >* hexBond 
    {new vector<vector<double> >(nbonds,vec)};  
    for (std::set<pair<int,int> >::iterator it = (*bonds)[ibin].begin();
	 it != (*bonds)[ibin].end(); it++) {
      ind1 = it->first;
      ind2 = it->second;
      xi = (*pos)[ibin][ind1][0];
      yi = (*pos)[ibin][ind1][1];
      xj = (*pos)[ibin][ind2][0];
      yj = (*pos)[ibin][ind2][1];
      dx = diff(lx,xi,xj);
      dy = diff(ly,yi,yj);
      (*posBond)[n][0] = xi-0.5*dx;
      (*posBond)[n][1] = yi-0.5*dy;
      theta = 6.0*(atan2(-dy,-dx)+pi);
      (*hexBond)[n][0] = cos(theta);
      (*hexBond)[n][1] = sin(theta);
      n++;
    }
    // Compute the correlation between bonds
    for (int i {}; i < nbonds; i++) {
      xi = (*posBond)[i][0];
      yi = (*posBond)[i][1];
      hexxi = (*hexBond)[i][0];
      hexyi = (*hexBond)[i][1];
      for (int j {}; j <= i; j++) {
	xj = (*posBond)[j][0];
	yj = (*posBond)[j][1];
	dr = dist(diff(lx,xi,xj),diff(ly,yi,yj));
	idbin = static_cast<int>(floor((dr-min)/binSize));
	hexxj = (*hexBond)[j][0];
	hexyj = (*hexBond)[j][1];
	corrx = hexxi*hexxj+hexyi*hexyj;
	corry = hexxj*hexyi-hexxi*hexyj;
	hexDistrb[idbin][0] += corrx;
	hexDistrb[idbin][1] += corry;
	distrbCount[idbin]++;
      }
    }
    delete posBond;
    delete hexBond;
  }

  // Output distribution
  double left, centre, right;
  for (int i {}; i < ndbins; i++) {
    // Normalise results
    if (distrbCount[i] > 0) {
      hexDistrb[i][0] /= static_cast<double>(distrbCount[i]);
      hexDistrb[i][1] /= static_cast<double>(distrbCount[i]);
    }
    left = i*binSize;
    right = (i+1)*binSize;
    centre = (left+right)*0.5;
    writer << left << " " << centre << " " << right << " " 
	   << hexDistrb[i][0] << " " << hexDistrb[i][1] << endl;
  }
  writer.close();

  delete pos;
  delete bonds;
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
