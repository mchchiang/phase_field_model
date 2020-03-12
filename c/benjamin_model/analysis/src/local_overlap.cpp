// local_overlap.cpp
// A program to compute the hexatic order parameter

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <omp.h>
#include "position.hpp"
#include "array.hpp"

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

int main(int argc, char* argv[]) {
  
  if (argc != 14) {
    cout << "usage: local_overlap npoints lx ly cellLx cellLy startTime "
	 << "endTime timeInc fieldTimeInc posFile fieldFileDir "
	 << "fieldFileName outFile"  << endl;
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
  string fieldFileDir {argv[++argi]};
  string fieldFileName {argv[++argi]};
  string outFile {argv[++argi]};
  
  PositionReader posReader;
  posReader.open(posFile, npoints, lx, ly, timeInc);
  if (!posReader.isOpen()) {
    cout << "Problem with opening position file!" << endl;
    return 1;
  }

  // Read the position data
  int nbins {static_cast<int>((endTime-startTime)/fieldTimeInc)+1};
  double*** pos {create3DDoubleArray(nbins, npoints, 2)};
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
	pos[ibin][i][0] = posReader.getPosition(i, 0);
	pos[ibin][i][1] = posReader.getPosition(i, 1);
      }
    }
  }
  posReader.close();  
  
  // Read cell phase fields
  ostringstream oss;
  ifstream fieldReader;
  istringstream iss;
  string line;
  double phi;
  double xavg, yavg, mass;
  double*** localFields {create3DDoubleArray(npoints, cellLx, cellLy)};
  int** pos0 {create2DIntArray(npoints, 2)};
  
  int i, j, m, n;
  int x, y, xn, yn;
  double olapAvg;
  double x0m, y0m, x0n, y0n, dx0mn, dy0mn;
  double phim, phin, area;
  vector<double> olap (npoints, 0.0);
  
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening the file: " << outFile << endl;
    return 1;
  }

  for (ibin = 0; ibin < nbins; ibin++) {
    time = ibin*fieldTimeInc+startTime;
    cout << "Doing bin " << ibin << endl;

#pragma omp parallel default(none) \
  shared(lx, ly, cellLx, cellLy, olap, localFields, pos, pos0, npoints, ibin, fieldFileDir, fieldFileName, time) \
  private(i, j, m, n, x, y, xavg, yavg, xn, yn, x0m, y0m, x0n, y0n, dx0mn, dy0mn, phi, phim, phin, area, mass, olapAvg, iss, oss, fieldReader, line)
    {
#pragma omp for schedule(static)
      for (n = 0; n < npoints; n++) {
	oss.str("");
	oss.clear();
	oss << fieldFileDir << "/cell_" << n << "/" 
	    << fieldFileName << "cell_" << n << ".dat." << time;
	fieldReader.open(oss.str());
	// Compute the local field centre of mass
	xavg = 0.0;
	yavg = 0.0;
	mass = 0.0;
	while (getline(fieldReader, line)) {
	  if (line.size() == 0) continue;
	  iss.clear();
	  iss.str(line);
	  iss >> x >> y >> phi;
	  localFields[n][x][y] = phi;
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
	x = wrap(lx, static_cast<int>(round(pos[ibin][n][0]-xavg)));
	y = wrap(ly, static_cast<int>(round(pos[ibin][n][1]-yavg)));
	pos0[n][0] = x;
	pos0[n][1] = y;
	fieldReader.close();
      }

    // Compute overlap
#pragma omp for schedule(static)
      for (m = 0; m < npoints; m++) {
	olapAvg = 0.0;
	area = 0.0;
	// Compute the area of the cell
	for (i = 0; i < cellLx; i++) {
	  for (j = 0; j < cellLy; j++) {
	    if (localFields[m][i][j] < 1.0) continue;
	    area += 1.0;
	  }
	}
	x0m = pos0[m][0];
	y0m = pos0[m][1];
	// Consider pairwise overlap
	for (n = 0; n < npoints; n++) {
	  if (m == n) continue; // Ignore self overlap
	  x0n = pos0[n][0];
	  y0n = pos0[n][1];
	  dx0mn = diff(lx, x0m, x0n);
	  dy0mn = diff(ly, y0m, y0n);
	  if (abs(dx0mn) > cellLx || abs(dy0mn) > cellLy) continue;
	  for (i = 0; i < cellLx; i++) {
	    xn = dx0mn+i;
	    if (xn < 0 || xn >= cellLx) continue;
	    for (j = 0; j < cellLy; j++) {
	      phim = localFields[m][i][j];
	      if (phim < 1.0) continue;
	      phim = 1.0; 
	      yn = dy0mn+j;
	      if (yn < 0 || yn >= cellLy) continue;
	      phin = localFields[n][xn][yn];
	      if (phin < 1.0) continue;
	      olapAvg += 1.0;
	    }
	  } // Close loop over local field of phi_n
	} // Close loop over all other cell n
	olapAvg /= area;
	olap[m] = olapAvg;
      } // Close loop over all cells
    } // Close parallel region
    
    // Output overlap data
    writer << "Cells: " << npoints << endl;
    writer << "Timestep: " << time << endl;
    writer << std::setprecision(5) << std::fixed;
    for (m = 0; m < npoints; m++) {
      writer << olap[m] << endl;
    }
    writer.unsetf(std::ios_base::floatfield);
  } // Close time loop
  writer.close();

  // Clean up resources
  free(pos);
  free(pos0);
  free(localFields);
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
