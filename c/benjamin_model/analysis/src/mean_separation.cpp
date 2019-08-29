// mean_separation.cpp
// A program that computes the mean separation between nearest neighbours

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <omp.h>
#include <cmath>
#include <cstdio>
#include "position.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::string;
using std::istringstream;

double getPeriodicDist(double len, double dr);

int main (int argc, char* argv[]) {
  
  if (argc != 12) {
    cout << "usage: mean_separation npoints lx ly startTime endTime "
	 << "timeInc endShiftTime posFile posBulkFile neighFile outFile" 
	 << endl;
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
  string neighFile {argv[++argi]};
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
  cout << "Reading position data ..." << endl;
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
  reader.close();
  cout << "Done reading position data" << endl;

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
  cout << "Reading bulk position data ..." << endl;
  for (int i = 0; i < nbins; i++) {
    cmReader >> t >> xcm >> ycm;
    totCM[i][0] = xcm;
    totCM[i][1] = ycm;
  }
  cmReader.close();
  cout << "Done reading bulk position data" << endl;

  // Read neighbour files
  vector<int> neighvec;
  neighvec.reserve(6); // Reserve space for at least 6 nearest neighbours
  vector<vector<int> > neighvecs (npoints, neighvec);
  vector<vector<vector<int> > >* neighbours
  {new vector<vector<vector<int> > >(nbins, neighvecs)};
  ifstream neighReader;
  neighReader.open(neighFile);
  if (!neighReader) {
    cout << "Error in reading neighbour file!" << endl;
    return 1;
  }
  cout << "Reading neighbour data ..." << endl;

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
  cout << "Done reading neighbour data" << endl;
  
  // Calculate mean separation
  int endShiftBin {static_cast<int>((endShiftTime-startTime)/timeInc)+1};
  if (endShiftBin >= nbins) endShiftBin = nbins-1;

  int i, j, k, l;
  int numOfNeighs;
  double x, y, dx, dy, localDistAvg, distAvg;
  //vector<double> meanSep (nbins, 0.0);
  //vector<double> count (nbins, 0.0);
  vector<long> count;
  vector<long> totCount (nbins, 0L);
  vector<double> meanSep;
  vector<double> totMeanSep (nbins, 0.0);

#pragma omp parallel default(none)					\
  shared(nbins, npoints, pos, totCM, neighbours, totMeanSep, totCount, endShiftBin, lx, ly) \
  private(ibin, i, j, k, l, numOfNeighs, x, y, dx, dy, localDistAvg, distAvg, meanSep, count)
{
  
  // Initialise arrays
  meanSep = vector<double>(nbins, 0.0);
  count = vector<long>(nbins, 0L);
    
#pragma omp for schedule(dynamic, 10)  
  for (i = 0; i < endShiftBin; i++) {
    for (j = i; j < nbins; j++) {
      ibin = abs(j-i);
      distAvg = 0.0;
      for (k = 0; k < npoints; k++) {
	x = (*pos)[i][k][0]-totCM[i][0];
	y = (*pos)[i][k][1]-totCM[i][1];
	localDistAvg = 0.0;
	numOfNeighs = static_cast<int>((*neighbours)[i][k].size());
	for (l = 0; l < numOfNeighs; l++) {
	  dx = x-((*pos)[j][(*neighbours)[i][k][l]][0]-totCM[j][0]);
	  dy = y-((*pos)[j][(*neighbours)[i][k][l]][1]-totCM[j][1]);
	  dx = getPeriodicDist(lx, dx);
	  dy = getPeriodicDist(ly, dy);
	  localDistAvg += sqrt(dx*dx+dy*dy);
	}
	distAvg += localDistAvg/static_cast<double>(numOfNeighs);
      }
      distAvg /= static_cast<double>(npoints);
      meanSep[ibin] += distAvg;
      count[ibin]++;
    }
  }
  printf("Done loop\n");
  for (i = 0; i < nbins; i++) {
#pragma omp atomic
    totMeanSep[i] += meanSep[i];
#pragma omp atomic
    totCount[i] += count[i];
  }
} // End of parallel region
 cout << "Done parallel computation" << endl;
  // Normalise
  for (i = 0; i < nbins; i++) {
    totMeanSep[i] /= static_cast<double>(totCount[i]);
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
	   << totMeanSep[i] << endl;
    writer.unsetf(std::ios_base::floatfield);
  }
  writer.close();
  
  delete pos;
  delete neighbours;
}


double getPeriodicDist(double len, double dr) {
  // Get the minimum separation
  dr = fabs(fmod(dr,len));
  double drAlter = fabs(len-dr);
  return drAlter > dr ? dr : drAlter;
}
