// asphericity.cpp
// A code to compute the asphericity of a 2D object using 
// the gyration tensor

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

double asphere(double gxx, double gyy, double gxy);

int main(int argc, char* argv[]) {
  
  if (argc != 6) {
    cout << "usage: asphericity npoints startTime endTime "
	 << "gyrFile outFile" << endl;
    return 1;
  }

  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  string gyrFile {argv[++argi]};
  string outFile {argv[++argi]};
  
  ifstream reader;
  reader.open(gyrFile);
  if (!reader) {
    cout << "Problem with opening gyration file!" << endl;
    return 1;
  }

  string line, str;
  istringstream iss;
  long time; 
  double gxx, gyy, gxy, b;
  double asphereAvg {};
  double asphereAvgSq {};
  long count {};
  while (!reader.eof()) {
    // Read the two header lines and get time
    getline(reader, line);
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
	iss >> gxx >> gyy >> gxy;
	b = asphere(gxx, gyy, gxy);
	asphereAvg += b;
	asphereAvgSq += b*b;
	count++;
      }
    }
  }
  reader.close();
  
  // Normalise
  double n {static_cast<double>(count)};
  asphereAvg /= n;
  asphereAvgSq /= n;
  double stdev {n/(n-1.0)*(asphereAvgSq-asphereAvg*asphereAvg)};
  double stderr {stdev/sqrt(n)};

  // Output results
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening output file!" << endl;
    return 1;
  }
  writer << std::setprecision(10) << std::fixed;
  writer << asphereAvg << " " << stdev << " " << stderr << endl;
  writer.close();
}

double asphere(double gxx, double gyy, double gxy) {
  double b {(gxx+gyy)/2.0};
  double c {gxx*gyy-gxy*gxy};
  double dis {b*b-c};
  double lam1 {b+sqrt(dis)};
  double lam2 {b-sqrt(dis)};
  return (lam1-lam2)/(lam1+lam2);
}
