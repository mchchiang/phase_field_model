// msd.cpp
// A code that computes the mean square displacement from the trajectory file

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include "position.hpp"

using std::cout;
using std::endl;
using std::ofstream;
using std::vector;
using std::map;
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
  reader.close();

  map<long, long long> count;
  map<long, double> r2avg;
  map<long, double> r4avg;
  vector<double> vec (2, 0.0);
  vector<vector<double> > dr (npoints, vec);
  vector<vector<double> > rold (npoints, vec);

  long time, dataTime;
  double x, y, dx, dy, r2, r4; 

  // Initial point is by definition zero
  r2avg[0] = 0.0;
  r4avg[0] = 0.0; 

  for (long t {startTime}; t < endTime; t += timeInc) {
    cout << "Doing t = " << t << endl;
    reader.open(posFile, npoints, lx, ly, timeInc);
    for (long i {startTime}; i < t; i += timeInc) {
      if (!reader.skipFrame()) break;
    }
    while (reader.nextFrame()) {
      time = reader.getTime();
      dataTime = time-t;
      if (time < t) {
	continue;
      } else if (time == t) {
	for (int i {}; i < npoints; i++) {
	  rold[i][0] = reader.getUnwrappedPosition(i, 0);
	  rold[i][1] = reader.getUnwrappedPosition(i, 1);
	  dr[i][0] = 0.0;
	  dr[i][1] = 0.0;
	} 
      } else if (time > t && time <= endTime) {
	for (int i {}; i < npoints; i++) {
	  x = reader.getUnwrappedPosition(i, 0);
	  y = reader.getUnwrappedPosition(i, 1);
	  dx = x - rold[i][0];
	  dy = y - rold[i][1];
	  dr[i][0] += dx;
	  dr[i][1] += dy;
	  r2 = (dr[i][0]*dr[i][0]+dr[i][1]*dr[i][1]);
	  r4 = r2*r2;
	  r2avg[dataTime] += r2;
	  r4avg[dataTime] += r4;
	  count[dataTime]++;
	  rold[i][0] = x;
	  rold[i][1] = y;
	}
      } else {
	break;
      }
    }
    reader.close();
  }

  // Normalise results
  for (map<long,double>::iterator it {r2avg.begin()}; 
       it != r2avg.end(); it++) {
    r2avg[it->first] /= static_cast<double>(count[it->first]);
    r4avg[it->first] /= static_cast<double>(count[it->first]);
  }

  // Output results
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening output file!" << endl;
    return 1;
  }
  writer << std::fixed << std::setprecision(5);

  for (map<long,double>::iterator it {r2avg.begin()}; 
       it != r2avg.end(); it++) {
    writer << it->first << " " << it->second << " " 
	   << r4avg[it->first] << endl;
  }
  writer.close();
}
