// defect_intermittence_signal.cpp
// A program to determine if the defect signal is intermittence or not

#include <iostream>
#include <fstream>
#include <sstream>
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
  
  if (argc != 8) {
    cout << "usage: defect_intermittence startTime endTime timeInc"
	 << "defectThres timeThres neighDiffFile outFile" << endl;
    return 1;
  }

  int argi {};
  long startTime {stol(string(argv[++argi]), nullptr, 10)};
  long endTime {stol(string(argv[++argi]), nullptr, 10)};
  long timeInc {stol(string(argv[++argi]), nullptr, 10)};
  double defectThres {stod(string(argv[++argi]), nullptr)};
  double timeThres {stod(string(argv[++argi]), nullptr)};
  string neighDiffFile {argv[++argi]};
  string outFile {argv[++argi]};
  
  // Read the defect signal
  ifstream reader;
  reader.open(neighDiffFile);
  if (!reader) {
    cout << "Problem with opening the defect signal file!" << endl;
    return 1;
  }

  istringstream iss;
  string line;
  long time;
  double defect;
  int nbins {static_cast<int>((endTime-startTime)/timeInc)+1};
  vector<int> state (nbins, 0);
  int ibin;
  
  // Convert the defect signal into a binary data set
  while (getline(reader, line)) {
    iss.clear();
    iss.str(line);
    iss >> time >> defect;
    if (time < startTime) {
      continue;
    } else if (time > endTime) {
      break;
    } else {
      ibin = static_cast<int>((time-startTime)/timeInc);
      if (defect > defectThres) {
	state[ibin] = 1;
      } else {
	state[ibin] = 0;
      }
    }
  }

  // Smooth the signal - remove any jumps within an interval less than 
  // a certain threshold (say 10^5 timesteps)
  int shift {static_cast<int>(timeThres/timeInc)};
  int currentState {state[0]};
  int prevState;
  for (int i {1}; i < nbins; i++) {
    prevState = currentState;
    currentState = state[i];
    if (currentState != prevState) {
      for (int j {i+1}; j <= (i+shift) && j < nbins; j++) {
	if (state[j] == prevState) {
	  state[i] = prevState;
	  currentState = prevState;
	  break;
	}
      }
    }
  }

  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  for (int i {}; i < nbins; i++) {
    time = startTime+timeInc*i;
    writer << time << " " << state[i] << endl;
  }
  writer.close();
}
