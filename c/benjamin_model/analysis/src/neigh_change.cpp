// neigh_change.cpp
// A program to check if there are neighbour changes

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <cmath>
#include "position.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::string;
using std::vector;
using std::set;

double dist(double x, double y);
double dot(double x1, double y1, double x2, double y2);
template <typename T> int sgn(T val);

int main(int argc, char* argv[]) {
  
  if (argc != 7) {
    cout << "usage: neighbour_change npoints startTime endTime "
	 << "timeInc neighFile outFile" << endl;
    return 1;
  }

  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  long timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  string neighFile {argv[++argi]};
  string outFile {argv[++argi]};
    
  // Read the neighbour data
  int nbins {static_cast<int>((endTime-startTime)/timeInc)+1};
  int ibin;
  long time;
  set<int> neighvec;
  vector<set<int> > neighvecs (npoints, neighvec);
  vector<vector<set<int> > >* neighbours
  {new vector<vector<set<int> > >(nbins, neighvecs)};
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
	    (*neighbours)[ibin][i].insert(neighIndex);
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
  
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  
  int count;
  int shift {1};
  for (ibin = shift; ibin < nbins; ibin++) {
    count = 0;
    for (int n {}; n < npoints; n++) {
      if ((*neighbours)[ibin][n] != (*neighbours)[ibin-shift][n]) {
	count++;
      }
    }
    writer << (ibin*timeInc+startTime) << " " << count << "\n";
  }
  writer.close();

  // Clean up resources
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
