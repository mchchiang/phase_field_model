// msd_time.cpp
// A code that computes the mean square displacement from the trajectory file

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include "position.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::vector;
using std::string;

int main (int argc, char* argv[]) {
  
  if (argc != 11) {
    cout << "usage: msd_time npoints lx ly startTime endTime "
	 << "timeInc timeShift posFile posBulkFile outFile" << endl;
    return 1;
  }
  
  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
  int lx {stoi(string(argv[++argi]), nullptr, 10)};
  int ly {stoi(string(argv[++argi]), nullptr, 10)};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  long timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  long timeShift {stoi(string(argv[++argi]), nullptr, 10)};
  string posFile {argv[++argi]};
  string posBulkFile {argv[++argi]};
  string outFile {argv[++argi]};
  
  PositionReader reader;
  reader.open(posFile, npoints, lx, ly, timeInc);
  if (!reader.isOpen()) {
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
  cout << "Reading data ..." << endl;
  while (reader.nextFrame()) {
    time = reader.getTime();
    if (time < startTime) {
      continue;
    } else if (time >= startTime && time <= endTime) {
      ibin = static_cast<int>((time-startTime)/timeInc);
      for (int i {}; i < npoints; i++) {
	(*pos)[ibin][i][0] = reader.getUnwrappedPosition(i, 0);
	(*pos)[ibin][i][1] = reader.getUnwrappedPosition(i, 1);
      }
    } else {
      break;
    }
  }
  cout << "Done reading data" << endl;
  reader.close();

  // Read bulk cm
  vector<vector<double> > totCM (nbins, vec);
  ifstream cmReader;
  cmReader.open(posBulkFile);
  if (!cmReader) {
    cout << "Error in reading bulk CM file!" << endl;
    return 1;
  }
  long t;
  double xcm, ycm;
  istringstream iss;
  string line;
  while (getline(cmReader, line)) {
    iss.clear();
    iss.str(line);
    iss >> t >> xcm >> ycm;
    if (t > endTime) {
      break;
    } else if (t < startTime) {
      continue;
    } else {
      ibin = static_cast<int>((t-startTime)/timeInc);
      totCM[ibin][0] = xcm;
      totCM[ibin][1] = ycm;
    }
  }
  cmReader.close();

  double dx, dy, r2; 
  double dxcm;
  double dycm;
  int shift {static_cast<int>(timeShift/timeInc)};

  // Output results
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening output file!" << endl;
    return 1;
  }
  
  for (ibin = shift; ibin < nbins; ibin++) {
    r2 = 0.0;
    dxcm = totCM[ibin][0]-totCM[ibin-shift][0];
    dycm = totCM[ibin][1]-totCM[ibin-shift][1];
    for (int n {}; n < npoints; n++) {
      dx = ((*pos)[ibin][n][0]-(*pos)[ibin-shift][n][0])-dxcm;
      dy = ((*pos)[ibin][n][1]-(*pos)[ibin-shift][n][1])-dycm;
      r2 = dx*dx + dy*dy;
    }
    r2 /= static_cast<double>(npoints);
    writer << (ibin*timeInc+startTime) << " " 
	   << std::setprecision(10) << std::fixed << r2 << "\n";
    writer.unsetf(std::ios_base::floatfield);    
  }
  writer.close();

  delete pos;
}
