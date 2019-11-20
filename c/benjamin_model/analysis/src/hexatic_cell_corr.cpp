// hexatic_cell_corr.cpp
// A program to compute the correlation of the hexatic order parameter

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
  
  if (argc != 13) {
    cout << "usage: hexatic_cell_corr npoints lx ly min max binSize startTime "
	 << "endTime timeInc posFile hexFile outFile" << endl;
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
  string hexFile {argv[++argi]};
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
  
  ifstream reader;
  reader.open(hexFile);
  if (!reader) {
    cout << "Problem with opening the hexatic file!" << endl;
    return 1;
  }

  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  writer << std::setprecision(5) << std::fixed;
  
  int idbin;
  int ndbins {static_cast<int>(ceil((max-min)/binSize))};
  vector<double> hexDistrb(ndbins, 0.0);
  vector<int> distrbCount(ndbins, 0);
  
  vector<vector<double> >* hexOrder 
  {new vector<vector<double> > (npoints, vec)};

  double hexx, hexy, mag, proj;
  double xi, xj, yi, yj, dr;
  double hexxi, hexxj, hexyi, hexyj, corrx, corry;
  string line, timestr;
  istringstream iss;
  while (getline(reader, line)) {
    // Read header lines
    if (!getline(reader, line)) break;
    iss.clear();
    iss.str(line);
    iss >> timestr >> time;
    if (time < startTime) {
      // Skip data before start time
      for (int i {}; i < npoints; i++) {
	if (!getline(reader, line)) break;
      }
    } else if (time > endTime) {
      break;
    } else {
      // Read the hexatic order for a single time frame
      ibin = static_cast<int>((time-startTime)/timeInc);
      for (int i {}; i < npoints; i++) {
	if (!getline(reader, line)) break;
	iss.clear();
	iss.str(line);
	iss >> hexx >> hexy >> mag >> proj;
	(*hexOrder)[i][0] = hexx;
	(*hexOrder)[i][1] = hexy;
      }
      // Compute the correlation
      for (int i {}; i < npoints; i++) {
	xi = (*pos)[ibin][i][0];
	yi = (*pos)[ibin][i][1];
	hexxi = (*hexOrder)[i][0];
	hexyi = (*hexOrder)[i][1];
	for (int j {}; j <= i; j++) {
	  xj = (*pos)[ibin][j][0];
	  yj = (*pos)[ibin][j][1];
	  dr = dist(diff(lx,xi,xj), diff(ly,yi,yj));
	  idbin = static_cast<int>(floor((dr-min)/binSize));
	  if (idbin < 0 || idbin > ndbins) {
	    cout << "Distance out of range: dr = " << dr << endl;
	    continue;
	  }
	  hexxj = (*hexOrder)[j][0];
	  hexyj = (*hexOrder)[j][1];
	  corrx = hexxi*hexxj+hexyi*hexyj;
	  corry = hexxj*hexyi-hexxi*hexyj;
	  hexDistrb[idbin] += dist(corrx,corry);
	  distrbCount[idbin]++;
	}
      }
    }
  }
  reader.close();

  // Output distribution
  double left, centre, right;
  for (int i {}; i < ndbins; i++) {
    // Normalise results
    if (distrbCount[i] > 0) {
      hexDistrb[i] /= static_cast<double>(distrbCount[i]);
    }
    left = i*binSize;
    right = (i+1)*binSize;
    centre = (left+right)*0.5;
    writer << left << " " << centre << " " << right << " " 
	   << hexDistrb[i] << " " 
	   << distrbCount[i]/(2.0*M_PI*left*binSize)/nbins << endl;
  }
  writer.close();

  delete pos;
  delete hexOrder;
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
