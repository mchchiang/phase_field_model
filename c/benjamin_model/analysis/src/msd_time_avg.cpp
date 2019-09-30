// msd.cpp
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
    cout << "usage: msd npoints lx ly startTime endTime "
	 << "timeInc endShiftTime posFile posBulkFile outFile" << endl;
    return 1;
  }
  
  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
  int lx {stoi(string(argv[++argi]), nullptr, 10)};
  int ly {stoi(string(argv[++argi]), nullptr, 10)};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  long timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  long endShiftTime {stoi(string(argv[++argi]), nullptr, 10)};
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

  int i, j, k;
  int endShiftBin {static_cast<int>((endShiftTime-startTime)/timeInc)+1};
  if (endShiftBin >= nbins) endShiftBin = nbins-1;
  double dx, dy, r2, r4; 
  vector<long> count;
  vector<long> totCount (nbins, 0);
  vector<double> r2avg, r4avg;
  vector<double> totR2Avg (nbins, 0.0);
  vector<double> totR4Avg (nbins, 0.0);
  double dxcm;
  double dycm;

#pragma omp parallel default(none) \
  shared(nbins, npoints, pos, totCM, totR2Avg, totR4Avg, totCount, endShiftBin) \
  private(time, ibin, i, j, k, dx, dy, r2, r4, r2avg, r4avg, count, dxcm, dycm)
{
    // Initialise array
    r2avg = vector<double>(nbins, 0.0);
    r4avg = vector<double>(nbins, 0.0);
    count = vector<long>(nbins, 0);

#pragma omp for schedule(dynamic,10)
    for (i = 0; i < endShiftBin; i++) {
      for (j = i+1; j < nbins; j++) {
	ibin = j-i;
	dxcm = totCM[j][0]-totCM[i][0];
	dycm = totCM[j][1]-totCM[i][1];
	for (k = 0; k < npoints; k++) {
	  dx = ((*pos)[j][k][0]-(*pos)[i][k][0])-dxcm;
	  dy = ((*pos)[j][k][1]-(*pos)[i][k][1])-dycm;
	  r2 = dx*dx + dy*dy;
	  r4 = r2*r2;
	  r2avg[ibin] += r2;
	  r4avg[ibin] += r4;
	}
	count[ibin] += npoints;
      }     
    }

    for (i = 0; i < nbins; i++) {
#pragma omp atomic
      totR2Avg[i] += r2avg[i];
#pragma omp atomic
      totR4Avg[i] += r4avg[i];
#pragma omp atomic
      totCount[i] += count[i];
    }
}
  // Normalise results
  for (i = 1; i < nbins; i++) { // first 
    totR2Avg[i] /= static_cast<double>(totCount[i]);
    totR4Avg[i] /= static_cast<double>(totCount[i]);
  }

  // Output results
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening output file!" << endl;
    return 1;
  }
  
  for (i = 0; i < nbins; i++) {
    writer << (i*timeInc) << " " << std::setprecision(10) << std::fixed 
	   << totR2Avg[i] << " " << totR4Avg[i] << endl;
    writer.unsetf(std::ios_base::floatfield);
  }
  writer.close();

  delete pos;
}
