// msd.cpp
// A code that computes the mean square displacement from the trajectory file

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "position.hpp"

using std::cout;
using std::endl;
using std::ofstream;
using std::vector;
using std::string;

int main (int argc, char* argv[]) {
  
  if (argc != 9) {
    cout << "usage: msd npoints lx ly startTime endTime "
	 << "timeInc posFile outFile" << endl;
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
  
  PositionReader reader;
  reader.open(posFile, npoints, lx, ly, timeInc);
  if (!reader.isOpen()) {
    cout << "Problem with opening position file!" << endl;
    return 1;
  }

  int nbins {static_cast<int>((endTime-startTime)/timeInc)+1};
  vector<double> r2avg (nbins, 0.0);
  vector<double> r4avg (nbins, 0.0);

  vector<double> vec (2, 0.0);
  vector<vector<double> > dr (npoints, vec);
  vector<vector<double> > rold (npoints, vec);
  
  // Do point average for now and don't do time average
  int ibin;
  long time;
  double x, y, dx, dy, r2, r4; 
  while (reader.nextFrame()) {
    time = reader.getTime();
    ibin = static_cast<int>((time-startTime) / timeInc);
    if (time < startTime) {
      continue;
    } else if (time == startTime) {
      for (int i {}; i < npoints; i++) {
	rold[i][0] = reader.getUnwrappedPosition(i, 0);
	rold[i][1] = reader.getUnwrappedPosition(i, 1);
      } 
    } else if (time > startTime && time <= endTime) {
      for (int i {}; i < npoints; i++) {
	x = reader.getUnwrappedPosition(i, 0);
	y = reader.getUnwrappedPosition(i, 1);
	dx = x - rold[i][0];
	dy = y - rold[i][1];
	dr[i][0] += dx;
	dr[i][1] += dy;
	r2 = (dr[i][0]*dr[i][0]+dr[i][1]*dr[i][1]);
	r4 = r2*r2;
	r2avg[ibin] += r2;
	r4avg[ibin] += r4;
	rold[i][0] = x;
	rold[i][1] = y;
      }
      r2avg[ibin] /= static_cast<double>(npoints);
      r4avg[ibin] /= static_cast<double>(npoints);
    } else {
      break;
    }
  }
  reader.close();

  // Output result
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening output file!" << endl;
    return 1;
  }
  writer << std::fixed << std::setprecision(5);

  for (int i {}; i < nbins; i++) {
    writer << (i*timeInc) << " " << r2avg[i] << " " << r4avg[i] << endl;
  }
  writer.close();
}
