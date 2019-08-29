// asphere_distrb.cpp
// A program to compute the distribution of asphericity values

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::string;
using std::vector;

int main(int argc, char* argv[]) {
  
  if (argc != 9) {
    cout << "usage: asphere_distrb npoints startTime endTime "
	 << "min max binSize shapeFile outFile" << endl;
    return 1;
  }

  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  double min {stod(string(argv[++argi]), nullptr)};
  double max {stod(string(argv[++argi]), nullptr)};
  double binSize {stod(string(argv[++argi]), nullptr)};
  string shapeFile {argv[++argi]};
  string outFile {argv[++argi]};
  
  ifstream reader;
  reader.open(shapeFile);
  if (!reader) {
    cout << "Problem with opening gyration file!" << endl;
    return 1;
  }
  
  // Set up distribution
  int nbins {static_cast<int>(ceil((max-min)/binSize))};
  vector<double> distrb (nbins, 0.0);
  double count {};

  string line, str;
  istringstream iss;
  long time; 
  int pixels;
  int ibin;
  double area, perimeter, pixelArea, chainPerimeter, asphere;
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
      // Compute asphericity
      for (int i {}; i < npoints; i++) {
	getline(reader, line);
	iss.clear();
	iss.str(line);
	iss >> perimeter >> area >> chainPerimeter >> pixelArea >> pixels;
	asphere = perimeter*perimeter/(4.0*M_PI*area);
	ibin = static_cast<int>((asphere-min)/binSize);
	if (ibin >= 0 && ibin < nbins) {
	  distrb[ibin] += 1.0;
	  count += 1.0;
	}
      }
    }
  }
  reader.close();

  // Normalise
  for (int i {}; i < nbins; i++) {
    distrb[i] /= count;
  }
  
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening output file!" << endl;
    return 1;
  }
  writer << std::setprecision(5) << std::fixed;
  
  // Output results to file
  double left, centre, right;
  for (int i {}; i < nbins; i++) {
    left = i*binSize+min;
    right = (i+1)*binSize+min;
    centre = (left+right)/2.0;
    writer << left << " " << centre << " " << right << " "
	   << distrb[i] << endl;	       
  }
  writer.close();
}
