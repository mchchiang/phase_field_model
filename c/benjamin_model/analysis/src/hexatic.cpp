// hexatic.cpp
// A program to compute the hexatic order parameter

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

double dist(double x, double y);
double dot(double x1, double y1, double x2, double y2);
template <typename T> int sgn(T val);

int main(int argc, char* argv[]) {
  
  if (argc != 11) {
    cout << "usage: hexatic npoints lx ly startTime endTime timeInc "
	 << "posFile neighFile cellOutFile avgOutFile" << endl;
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
  string neighFile {argv[++argi]};
  string cellOutFile {argv[++argi]};
  string avgOutFile {argv[++argi]};
  
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
  
  // Read the neighbour data
  vector<int> neighvec;
  neighvec.reserve(6); // Reserve space for at least 6 nearest neighbours
  vector<vector<int> > neighvecs (npoints, neighvec);
  vector<vector<vector<int> > >* neighbours
  {new vector<vector<vector<int> > >(nbins, neighvecs)};
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
	    (*neighbours)[ibin][i].push_back(neighIndex);
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
  
  ofstream cellWriter;
  cellWriter.open(cellOutFile);
  if (!cellWriter) {
    cout << "Problem with opening individual cell output file!" << endl;
    return 1;
  }
  ofstream avgWriter;
  avgWriter.open(avgOutFile);
  if (!avgWriter) {
    cout << "Problem with opening average output file!" << endl;
    return 1;
  }

  // Compute the local hexatic order for each cell at each time frame
  // phi_6i = 1/(N_i)\sum_{j=1}^{N_i} exp(i6\theta_{ij})
  vector<vector<double> > hexOrder (npoints, vec);
  vector<double> avgHexOrder (2, 0.0);
  double magHexOrder, avgMagHexOrder;
  double theta;
  double pi {M_PI};
  double x, y;
  double dx, dy, dx1, dx2, dy1, dy2;
  int numOfNeighs, ineigh;
  for (ibin = 0; ibin < nbins; ibin++) {
    time = ibin*timeInc;
    avgHexOrder[0] = 0.0;
    avgHexOrder[1] = 0.0;
    avgMagHexOrder = 0.0;
    for (int i {}; i < npoints; i++) {
      hexOrder[i][0] = 0.0;
      hexOrder[i][1] = 0.0;
      numOfNeighs = (*neighbours)[ibin][i].size();
      x = (*pos)[ibin][i][0];
      y = (*pos)[ibin][i][1];
      for (int j {}; j < numOfNeighs; j++) {
	ineigh = (*neighbours)[ibin][i][j];
	dx1 = x-(*pos)[ibin][ineigh][0];
	dy1 = y-(*pos)[ibin][ineigh][1];
	dx2 = lx-fabs(dx1);
	dy2 = ly-fabs(dy1);
	dx = (dx2 > fabs(dx1)) ? dx1 : -dx2*sgn(dx1);
	dy = (dy2 > fabs(dy1)) ? dy1 : -dy2*sgn(dy1);
	theta = atan2(-dy,-dx)+pi;
	hexOrder[i][0] += cos(6.0*theta);
	hexOrder[i][1] += sin(6.0*theta);
      }
      hexOrder[i][0] /= numOfNeighs;
      hexOrder[i][1] /= numOfNeighs;
      magHexOrder = dist(hexOrder[i][0], hexOrder[i][1]);
      avgMagHexOrder += magHexOrder;
      avgHexOrder[0] += hexOrder[i][0];
      avgHexOrder[1] += hexOrder[i][1];  
    }
    avgHexOrder[0] /= npoints;
    avgHexOrder[1] /= npoints;
    avgMagHexOrder /= npoints; 
    
    double magAvgHexOrder {dist(avgHexOrder[0], avgHexOrder[1])};

    // Write individual cell data
    cellWriter << "Cells: " << npoints << endl;
    cellWriter << "Timestep: " << time << endl;
    cellWriter << std::setprecision(5) << std::fixed;
    for (int i {}; i < npoints; i++) {
      // Compute projection to sample/global hexatic order
      magHexOrder = dist(hexOrder[i][0], hexOrder[i][1]);
      double projection {dot(hexOrder[i][0], hexOrder[i][1], avgHexOrder[0], 
			     avgHexOrder[1])/ (magHexOrder*magAvgHexOrder)};
      cellWriter << hexOrder[i][0] << " " << hexOrder[i][1] << " " 
		 << magHexOrder << " " << projection << endl;
    }
    cellWriter.unsetf(std::ios_base::floatfield);
    
    // Write average data
    avgWriter << time << " " << std::setprecision(10) << std::fixed
	      << avgHexOrder[0] << " " << avgHexOrder[1] << " " 
	      << magAvgHexOrder << " " << avgMagHexOrder << endl;
    avgWriter.unsetf(std::ios_base::floatfield);
  }
  cellWriter.close();
  avgWriter.close();

  // Clean up resources
  delete pos;
  delete neighbours;
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
