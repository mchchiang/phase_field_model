// self_int_scatter.cpp
// A code that computes the self-intermediate scattering function

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "position.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::string;

int main (int argc, char* argv[]) {
  
  if (argc != 13) {
    cout << "usage: self_int_scatter npoints lx ly nqvec r0 startTime endTime "
	 << "timeInc endShiftTime posFile posBulkFile outFile" << endl;
    return 1;
  }
  
  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
  int lx {stoi(string(argv[++argi]), nullptr, 10)};
  int ly {stoi(string(argv[++argi]), nullptr, 10)};
  int nqvec {stoi(string(argv[++argi]), nullptr, 10)};
  double r0 {stod(string(argv[++argi]), nullptr)};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  long timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  long endShiftTime {stoi(string(argv[++argi]), nullptr, 10)};
  string posFile {argv[++argi]};
  string posBulkFile {argv[++argi]};
  string outFile {argv[++argi]};

  double qmag {M_PI/r0};
  
  PositionReader reader;
  reader.open(posFile, npoints, lx, ly, timeInc);
  if (!reader.isOpen()) {
    cout << "Problem with opening position file!" << endl;
    return 1;
  }

  // Read the position data
  int nbins {static_cast<int>((endTime-startTime)/timeInc)+1};
  vector<double> vec (2, 0.0);
  vector<vector<double> > vec2 (nbins, vec);
  vector<vector<vector<double> > >* pos 
  {new vector<vector<vector<double > > >(npoints, vec2)};
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
	(*pos)[i][ibin][0] = reader.getUnwrappedPosition(i, 0);
	(*pos)[i][ibin][1] = reader.getUnwrappedPosition(i, 1);
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
  for (int i {}; i < nbins; i++) {
    cmReader >> t >> xcm >> ycm;
    totCM[i][0] = xcm;
    totCM[i][1] = ycm;
  }

  // Store a set of orientations for k to average over
  vector<vector<double> > qvec (nqvec, vec);
  double theta;
  for (int i {}; i < nqvec; i++) {
    theta = i*(2.0*M_PI/nqvec);
    qvec[i][0] = qmag*cos(theta);
    qvec[i][1] = qmag*sin(theta);
  }

  // Compute the self intermediate scattering function
  int i, j, n, m;
  vector<long> count;
  vector<long> totCount (nbins, 0L);
  vector<double> cosAvg, sinAvg;
  vector<double> totCosAvg (nbins, 0.0);
  vector<double> totSinAvg (nbins, 0.0);
  double dxcm, dycm, qdr;
  vector<vector<double> > dr;
  
  int endShiftBin {static_cast<int>((endShiftTime-startTime)/timeInc)+1};
  if (endShiftBin >= nbins) endShiftBin = nbins-1;
  
#pragma omp parallel default(none) \
  shared(nbins, npoints, pos, totCM, vec, nqvec, qvec, totCosAvg, totSinAvg) \
  shared(totCount, endShiftBin) \
  private(time, ibin, i, j, n, m, qdr, dxcm, dycm, dr, count, cosAvg, sinAvg)
{

  // Initialise the arrays
  cosAvg = vector<double>(nbins, 0.0);
  sinAvg = vector<double>(nbins, 0.0);
  count = vector<long>(nbins, 0L);
  dr = vector<vector<double> >(npoints, vec);

#pragma omp for schedule(dynamic, 10)
  for (n = 0; n < endShiftBin; n++) {
    for (m = n; m < nbins; m++) {
      ibin = m-n;
      dxcm = totCM[n][0]-totCM[m][0];
      dycm = totCM[n][1]-totCM[m][1];
      for (i = 0; i < npoints; i++) {
	dr[i][0] = (*pos)[i][n][0]-(*pos)[i][m][0]-dxcm;
	dr[i][1] = (*pos)[i][n][1]-(*pos)[i][m][1]-dycm;
      }
      for (j = 0; j < nqvec; j++) {
	for (i = 0; i < npoints; i++) {
	  qdr = qvec[j][0]*dr[i][0]+qvec[j][1]*dr[i][1];
	  cosAvg[ibin] += cos(qdr);
	  sinAvg[ibin] += sin(qdr);
	  count[ibin]++;
	}
      }
    }
  }
  
  // Add all the reuslts together
  for (i = 0; i < nbins; i++) {
#pragma omp atomic
    totCosAvg[i] += cosAvg[i];
#pragma omp atomic
    totSinAvg[i] += sinAvg[i];
#pragma omp atomic
    totCount[i] += count[i];
  }
}
  
  // Normalise
  for (i = 0; i < nbins; i++) {
    totCosAvg[i] /= static_cast<double>(totCount[i]);
    totSinAvg[i] /= static_cast<double>(totCount[i]);
  } 

  // Output results
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening output file!" << endl;
    return 1;
  }
  writer << std::fixed << std::setprecision(10);

  for (i = 0; i < nbins; i++) {
    writer << (i*timeInc) << " " << totCosAvg[i] << " " 
	   << totSinAvg[i] << endl;
  }
  writer.close();

  delete pos;
}
