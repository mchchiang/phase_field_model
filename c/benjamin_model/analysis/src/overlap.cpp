// overlap.cpp
// A program to compute the overlap parameter

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::string;

int main(int argc, char* argv[]) {
  
  if (argc != 8) {
    cout << "usage: overlap npoints lx ly startTime endTime "
	 << "shapeFile outFile" << endl;
    return 1;
  }

  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
  int lx {stoi(string(argv[++argi]), nullptr, 10)};
  int ly {stoi(string(argv[++argi]), nullptr, 10)};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  string shapeFile {argv[++argi]};
  string outFile {argv[++argi]};
  
  ifstream reader;
  reader.open(shapeFile);
  if (!reader) {
    cout << "Problem with opening gyration file!" << endl;
    return 1;
  }

  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening output file!" << endl;
    return 1;
  }
  
  string line, str;
  istringstream iss;
  long time; 
  int pixels;
  double area, perimeter, pixelArea, chainPerimeter;
  double totalCellArea {};
    while (getline(reader, line)) {
    // Read the two header lines and get time
    getline(reader, line);
    iss.clear();
    iss.str(line);
    iss >> str >> time;
    
    // Only use the data from the specified time period
    if (time < startTime) {
      // Skip the data in that time frame
      for (int i {}; i < npoints; i++) {
	getline(reader, line);
      }
    } else if (time > endTime) {
      break;
    } else {
      // Compute overlap
      totalCellArea = 0.0;
      for (int i {}; i < npoints; i++) {
	getline(reader, line);
	iss.clear();
	iss.str(line);
	iss >> perimeter >> area >> chainPerimeter >> pixelArea >> pixels;
	totalCellArea += pixels;
      }

      // Output results to file
      writer << time << " " << std::setprecision(10) << std::fixed
	     << (totalCellArea/(lx*ly)-1.0) << endl;
      writer.unsetf(std::ios_base::floatfield);
    }
  }
  reader.close();
  writer.close();
}
