// neigh_diff.cpp
// A program to compute the average difference from perfect 6-fold coordination

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

int main(int argc, char* argv[]) {
  
  if (argc != 6) {
    cout << "usage: neigh_diff npoints startTime endTime "
	 << "neighFile outFile" << endl;
    return 1;
  }

  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  string neighFile {argv[++argi]};
  string outFile {argv[++argi]};
    
  // Read the neighbour data
  ifstream neighReader;
  neighReader.open(neighFile);
  if (!neighReader) {
    cout << "Problem with opening neighbour file!" << endl;
    return 1;
  }

  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }

  string line, timestr;
  long time;
  istringstream iss;
  int count, totalCount, neighIndex;
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
	totalCount = 0;
	for (int i {}; i < npoints; i++) {
	  count = 0;
	  getline(neighReader, line);
	  iss.clear();
	  iss.str(line);
	  while (iss) {
	    iss >> neighIndex;
	    if (!iss) break;
	    count++;
	  }
	  totalCount += abs(count-6);
	}
	writer << time << " " << totalCount/static_cast<double>(npoints)
	       << endl;
      } else {
	break;
      }
    } else {
      break;
    }
  }
  neighReader.close();
  writer.close();
}
