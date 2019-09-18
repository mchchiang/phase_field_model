// colour_field.cpp
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

int main(int argc, char* argv[]) {
  
  if (argc != 9) {
    cout << "usage: colour_field npoints dataCol startTime endTime timeInc "
	 << "dataFile indexFieldFile outFieldFile" << endl;
    return 1;
  }

  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
  int dataCol {stoi(string(argv[++argi]), nullptr, 10)};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  long timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  string dataFile {argv[++argi]};
  string indexFieldFile {argv[++argi]};
  string outFieldFile {argv[++argi]};
  
  // Read data file
  ifstream dataReader;
  dataReader.open(dataFile);
  if (!dataReader) {
    cout << "Problem with opening the data file!" << endl;
    return 1;
  }
  
  int x, y, index;
  double value;
  string str, line, filename;
  istringstream iss;
  ostringstream oss;
  long time;
  vector<double> data (npoints, 0.0);
  ifstream indexFieldReader;
  ofstream outFieldWriter;
  while (getline(dataReader, line)) {
    if (!getline(dataReader, line)) break;
    // Get time
    iss.clear();
    iss.str(line);
    iss >> str >> time;
    if (time < startTime || (time-startTime)%timeInc != 0) {
      // Skip current time frame
      for (int i {}; i < npoints; i++)
	getline(dataReader, line);
    } else if (time > endTime) {
      break;
    } else {
      for (int i {}; i < npoints; i++) {
	getline(dataReader, line);
	iss.clear();
	iss.str(line);
	bool foundValue = true;
	for (int j {}; j <= dataCol; j++) {
	  if (!iss) {
	    foundValue = false;
	    break;
	  }
	  iss >> value;
	}
	if (foundValue) {
	  data[i] = value;
	}
      }
      // Read index field file
      oss.clear();
      oss.str("");
      oss << indexFieldFile << "." << time;
      indexFieldReader.open(oss.str());
      if (!indexFieldReader) continue;
      oss.clear();
      oss.str("");
      oss << outFieldFile << "." << time;
      outFieldWriter.open(oss.str());
      if (!outFieldWriter) continue;
      while (getline(indexFieldReader, line)) {
	if (line.size() > 0) {
	  iss.clear();
	  iss.str(line);
	  iss >> x >> y >> index;
	  outFieldWriter << x << " " << y << " " << std::setprecision(10) 
			 << std::fixed  << data[index-1] << endl;
	  outFieldWriter.unsetf(std::ios_base::floatfield);
	} else {
	  outFieldWriter << endl;
	}
      }
      indexFieldReader.close();
      outFieldWriter.close();
    }
  }
}
