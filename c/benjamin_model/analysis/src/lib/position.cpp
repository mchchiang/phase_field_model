// PositionReader.cpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "position.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::vector;
using std::string;

// Destructor
PositionReader::~PositionReader() {
  deleteData();
}

void PositionReader::initData() {
  vector<double> zeroDouble (2, 0.0);
  vector<int> zeroInt (2, 0);
  position = new vector<vector<double> >(npoints, zeroDouble);
  boundaryCount = new vector<vector<int> >(npoints, zeroInt);
}

void PositionReader::deleteData() {
  if (!position) delete position;
  if (!boundaryCount) delete boundaryCount;
}

bool PositionReader::open(string posFile, int npts,
			  double lx, double ly, int tInc) {
  if (fileOpen) {
    cout << "Close current file before opening another file" << endl;
    return false;
  }
  reader.open(posFile);
  if (!reader) {
    cout << "Problem with opening position file!" << endl;
    return false;
  }
  fileOpen = true;
  npoints = npts;
  boxSize = {lx, ly};
  timeInc = tInc;
  time = 0;
  initData();
  return true;
}

void PositionReader::close() {
  if (reader.is_open()) {
    reader.close();
  }
  deleteData();
  fileOpen = false;
}

double PositionReader::getPosition(int index, int comp) const {
  return (*position)[index][comp];
}

double PositionReader::getUnwrappedPosition(int index, int comp) const {
  return (*position)[index][comp] + 
    boxSize[comp] * (*boundaryCount)[index][comp];
}

bool PositionReader::nextFrame() {
  if (!fileOpen || reader.eof()) {
    return false;
  }
  string line, strTime;
  long newTime;
  double x, y;
  int ix, iy;
  istringstream iss;
  
  // Ignore header lines
  if (!getline(reader, line)) return false;
  if (!getline(reader, line)) return false;
  iss.clear();
  iss.str(line);
  iss >> strTime >> newTime;

  // Read bead position data
  for (int i {}; i < npoints; i++) {
    if (!getline(reader, line)) return false;
    iss.clear();
    iss.str(line);
    iss >> x >> y >> ix >> iy;
    (*position)[i][0] = x;
    (*position)[i][1] = y;
    (*boundaryCount)[i][0] = ix;
    (*boundaryCount)[i][1] = iy;
  }

  time = newTime;
  return true;
}

bool PositionReader::skipFrame() {
  if (!fileOpen || reader.eof()) {
    return false;
  }
  
  string line, strTime;
  long newTime;
  istringstream iss;
  // Ignore header lines
  if (!getline(reader, line)) return false;
  if (!getline(reader, line)) return false;
  iss.clear();
  iss.str(line);
  iss >> strTime >> newTime;
  for (int i {}; i < npoints; i++) {
    if (!getline(reader, line)) return false;
  }
  time = newTime;
  return true;
}

long PositionReader::getTime() const {
  return time;
}

bool PositionReader::isOpen() {
  return fileOpen;
}
