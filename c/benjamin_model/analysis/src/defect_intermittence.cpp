// defect_intermittence.cpp
// A program to determine if the defect signal is intermittence or not

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::string;

int main(int argc, char* argv[]) {
  
  if (argc != 7) {
    cout << "usage: defect_intermittence startTime endTime "
	 << "defectThres jumpThres timeFracThres neighDiffFile" << endl;
    return 1;
  }

  int argi {};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  double defectThres {stod(string(argv[++argi]), nullptr)};
  int jumpThres {stoi(string(argv[++argi]), nullptr, 10)};
  double timeFracThres {stod(string(argv[++argi]), nullptr)};
  string neighDiffFile {argv[++argi]};
  
  // Read the defect signal
  ifstream reader;
  reader.open(neighDiffFile);
  if (!reader) {
    cout << "Problem with opening the defect signal file!" << endl;
    return 1;
  }

  istringstream iss;
  string line;
  long time, prevTime {startTime}, solidTime {}, fluidTime {};
  double defect;
  int njumps {};
  int prevState {}, state {}; // 0 = solid, 1 = fluid
  
  while (getline(reader, line)) {
    iss.clear();
    iss.str(line);
    iss >> time >> defect;
    if (time < startTime) {
      continue;
    } else if (time > endTime) {
      break;
    } else {
      if (defect > defectThres) {
	state = 1;
      } else {
	state = 0;
      }
      if (time > startTime) {
	if (state == 1) {
	  fluidTime += (time-prevTime);
	} else {
	  solidTime += (time-prevTime);
	}
	if (state != prevState) {
	  njumps++;
	}
      }
      prevState = state;
    }
    prevTime = time;
  }

  double totalTime {static_cast<double>(endTime-startTime)};
  double fluidFrac {fluidTime/totalTime};
  double solidFrac {solidTime/totalTime};

  if (njumps > jumpThres && fluidFrac > timeFracThres && 
      solidFrac > timeFracThres) {
    cout << 1 << endl;
  } else {
    if (fluidFrac > solidFrac) {
      cout << 2 << endl;
    } else {
      cout << 0 << endl;
    }
  }
}
