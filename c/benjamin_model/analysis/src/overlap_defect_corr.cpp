// overlap_defect_corr.cpp
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
using std::ostringstream;
using std::string;
using std::vector;

double dist(double x, double y);
double dot(double x1, double y1, double x2, double y2);
template <typename T> int sgn(T val);
int wrap(int lx, int x);
int diff(int lx, int x1, int x2);

struct Defect {
  double x;
  double y;
  int cellIndex;
  int sign;
};

int main(int argc, char* argv[]) {
  
  if (argc != 14) {
    cout << "usage: overlap_defect_corr npoints lx ly cellLx cellLy startTime "
	 << "endTime timeInc fieldTimeInc posFile neighFile fieldFileRoot "
	 << "outFile" 
	 << endl;
    return 1;
  }

  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
  int lx {stoi(string(argv[++argi]), nullptr, 10)};
  int ly {stoi(string(argv[++argi]), nullptr, 10)};
  int cellLx {stoi(string(argv[++argi]), nullptr, 10)};
  int cellLy {stoi(string(argv[++argi]), nullptr, 10)};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  long timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  long fieldTimeInc {stoi(string(argv[++argi]), nullptr, 10)};
  string posFile {argv[++argi]};
  string neighFile {argv[++argi]};
  string fieldFileRoot {argv[++argi]};
  string outFile {argv[++argi]};
  
  PositionReader posReader;
  posReader.open(posFile, npoints, lx, ly, timeInc);
  if (!posReader.isOpen()) {
    cout << "Problem with opening position file!" << endl;
    return 1;
  }

  // Read the position data
  int nbins {static_cast<int>((endTime-startTime)/fieldTimeInc)+1};
  vector<double> vec (2, 0.0);
  vector<vector<double> > vec2 (npoints, vec);
  vector<vector<vector<double> > >* pos 
  {new vector<vector<vector<double > > >(nbins, vec2)};
  long time;
  int ibin;
  while (posReader.nextFrame()) {
    time = posReader.getTime();
    if (time > endTime) {
      break;
    } else if (time < startTime || (time-startTime)%fieldTimeInc != 0) {
      continue;
    } else {
      ibin = static_cast<int>((time-startTime)/fieldTimeInc);
      for (int i {}; i < npoints; i++) {
	// Use wrapped position here as we only care about relative distances
	// and orientations
	(*pos)[ibin][i][0] = posReader.getPosition(i, 0);
	(*pos)[ibin][i][1] = posReader.getPosition(i, 1);
      }
    }
  }
  posReader.close();
  
  // Read the neighbour data
  ifstream neighReader;
  neighReader.open(neighFile);
  if (!neighReader) {
    cout << "Problem with opening neighbour file!" << endl;
    return 1;
  }
  string line, timestr;
  time = 0L;
  vector<vector<Defect> >* defects {new vector<vector<Defect> >(nbins)};
  istringstream iss;
  int neighIndex, count;
  while (getline(neighReader, line)) {
    // Read time
    if (getline(neighReader, line)) {
      iss.clear();
      iss.str(line);
      iss >> timestr >> time;
      if (time > endTime) {
	break;
      } else if (time < startTime || (time-startTime)%fieldTimeInc != 0) {
	// Skip over unused data
	for (int n {}; n < npoints; n++) {
	  if (!getline(neighReader, line)) break;
	}
      } else {
	ibin = static_cast<int>((time-startTime)/fieldTimeInc);
	for (int n {}; n < npoints; n++) {
	  getline(neighReader, line);
	  iss.clear();
	  iss.str(line);
	  count = 0;
	  while (iss) {
	    iss >> neighIndex;
	    if (!iss) break;
	    count++;
	  }
	  if ((count-6) != 0) {
	    Defect d {(*pos)[ibin][n][0],(*pos)[ibin][n][1], n, count-6};
	    (*defects)[ibin].push_back(d);
	  }
	}
      }
    } else {
      break;
    }
  }
  neighReader.close();
  
  // Read cell phase fields
  ostringstream oss;
  ifstream fieldReader;
  double x, y, phi;
  double xavg, yavg, mass;
  vector<double> field1D (cellLy, 0.0);
  vector<vector<double> > field2D (cellLx, field1D);
  vector<vector<vector<double> > >* localFields
  {new vector<vector<vector<double> > >(npoints, field2D)};
  vector<vector<double> >* pos0 {new vector<vector<double> >(npoints, vec)};
  double positiveDefectOlapAvg {};
  double negativeDefectOlapAvg {};
  double totalOlapAvg {};
  double numOfPoints {};
  double numOfPositiveDefects {};
  double numOfNegativeDefects {};
  for (ibin = 0; ibin < nbins; ibin++) {
    // Skip time frames with no defects
    time = ibin*fieldTimeInc+startTime;
    if ((*defects)[ibin].size() == 0) continue;
    for (int n {}; n < npoints; n++) {
      oss.str("");
      oss.clear();
      oss << fieldFileRoot << "cell_" << n << ".dat." << time;
      fieldReader.open(oss.str());
      if (!fieldReader) {
	cout << "Problem with opening the file: " << oss.str() << endl;
	return 1;
      }
      // Compute the local field centre of mass
      xavg = 0.0;
      yavg = 0.0;
      mass = 0.0;
      while (getline(fieldReader, line)) {
	if (line.size() == 0) continue;
	iss.clear();
	iss.str(line);
	iss >> x >> y >> phi;
	(*localFields)[n][x][y] = phi;
	xavg += (phi*(x+0.5));
	yavg += (phi*(y+0.5));
	mass += phi;
      }
      if (mass > 0.0) {
	xavg /= mass;
	yavg /= mass;
      } else {
	xavg = 0.0;
	yavg = 0.0;
      }
      // Determine top left coordinates of the cell
      x = wrap(lx, static_cast<int>(round((*pos)[ibin][n][0]-xavg)));
      y = wrap(ly, static_cast<int>(round((*pos)[ibin][n][1]-yavg)));
      (*pos0)[n][0] = x;
      (*pos0)[n][1] = y;
      fieldReader.close();
    }
    
    // Compute overlap
    int xn, yn;
    double olapAvg;
    double x0m, y0m, x0n, y0n, dx0mn, dy0mn;
    double phim, phin, prod, area;
    vector<double> olap (npoints, 0.0);
    for (int m {}; m < npoints; m++) {
      olapAvg = 0.0;
      area = 0.0;
      // Compute the area of the cell
      for (int i {}; i < cellLx; i++) {
	for (int j {}; j < cellLy; j++) {
	  if ((*localFields)[m][i][j] < 1.0) continue;
	  area += 1.0;
	}
      }
      x0m = (*pos0)[m][0];
      y0m = (*pos0)[m][1];
      // Consider pairwise overlap
      for (int n {}; n < npoints; n++) {
	if (m == n) continue; // Ignore self overlap
	x0n = (*pos0)[n][0];
	y0n = (*pos0)[n][1];
	dx0mn = diff(lx, x0m, x0n);
	dy0mn = diff(ly, y0m, y0n);
	if (abs(dx0mn) > lx || abs(dy0mn) > ly) continue;
	for (int i {}; i < cellLx; i++) {
	  xn = dx0mn+i;
	  if (xn < 0 || xn >= cellLx) continue;
	  for (int j {}; j < cellLy; j++) {
	    phim = (*localFields)[m][i][j];
	    if (phim < 1.0) continue;
	    phim = 1.0; 
	    yn = dy0mn+j;
	    if (yn < 0 || yn >= cellLy) continue;
	    phin = (*localFields)[n][xn][yn];
	    phin = phin < 1.0 ? 0.0 : 1.0;
	    prod = phim*phin;
	    olapAvg += prod;
	  }
	} // Close loop over local field of phi_m
      } // Close loop over all other cell n
      olapAvg /= area;
      olap[m] = olapAvg;
      numOfPoints += 1.0;
      totalOlapAvg += olapAvg;
    } // Close loop over all cells
    for (auto& d : (*defects)[ibin]) {
      if (d.sign > 0) {
	cout << time << " +ve " << d.cellIndex << endl;
	numOfPositiveDefects++;
	positiveDefectOlapAvg += olap[d.cellIndex];
      } else {
	cout << time << " -ve " << d.cellIndex << endl;
	numOfNegativeDefects++;
	negativeDefectOlapAvg += olap[d.cellIndex];
      }
    }
  } // Close time loop

  // Normalise results
  totalOlapAvg /= static_cast<double>(numOfPoints);
  positiveDefectOlapAvg /= static_cast<double>(numOfPositiveDefects);
  negativeDefectOlapAvg /= static_cast<double>(numOfNegativeDefects);

  cout << totalOlapAvg << " " << positiveDefectOlapAvg << " " 
       << negativeDefectOlapAvg << endl;

  // Clean up resources
  delete pos;
  delete pos0;
  delete defects;
  delete localFields;
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

int wrap(int lx, int x) {
  int remainder {x%lx};
  if (remainder >= 0) return remainder;
  return lx + remainder;
}

int diff(int lx, int x1, int x2) {
  double dx1 = x1-x2;
  double dx2 = -sgn(dx1)*(lx-abs(dx1));
  if (abs(dx1) < abs(dx2)) return dx1;
  return dx2;
}
