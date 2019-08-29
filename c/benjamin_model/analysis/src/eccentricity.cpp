// eccentricity.cpp
// A program to compute the average eccentricity using the gyration tensor

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

double eccent(double gxx, double gyy, double gxy);

int main(int argc, char* argv[]) {
  
  if (argc != 6) {
    cout << "usage: eccentricity npoints startTime endTime "
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

  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening output file!" << endl;
    return 1;
  }

  string line, str;
  istringstream iss;
  long time; 
  double gxx, gyy, gxy, b;
  double eccentAvg {};
  double eccentAvgSq {};
  double n = static_cast<double>(npoints);
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
      // Compute eccentricity
      eccentAvg = 0.0;
      eccentAvgSq = 0.0;
      for (int i {}; i < npoints; i++) {
	getline(reader, line);
	iss.clear();
	iss.str(line);
	iss >> gxx >> gyy >> gxy;
	b = eccent(gxx, gyy, gxy);
	eccentAvg += b;
	eccentAvgSq += b*b;
      }
      
      // Normalise
      eccentAvg /= n;
      eccentAvgSq /= n;
      double stdev {n/(n-1.0)*(eccentAvgSq-eccentAvg*eccentAvg)};
      double stderr {stdev/sqrt(n)};
      
      // Output results
      writer << time << " " << std::setprecision(10) << std::fixed 
	     << eccentAvg << " " << stdev << " " << stderr << endl;
      writer.unsetf(std::ios_base::floatfield);
    }
  }
  reader.close();
  writer.close();
}

double eccent(double gxx, double gyy, double gxy) {
  double b {(gxx+gyy)/2.0};
  double c {gxx*gyy-gxy*gxy};
  double dis {b*b-c};
  double lam1 {b+sqrt(dis)};
  double lam2 {b-sqrt(dis)};
  return (lam1-lam2)/(lam1+lam2);
}
