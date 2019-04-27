/* position.hpp
 * A code that reads dump/position files.
 */

#ifndef POSITION_HPP
#define POSITION_HPP

#include <iostream>
#include <vector>
#include <string>

using std::vector;
using std::string;
using std::ifstream;

class PositionReader {

private:
  
  int npoints {};
  vector<double> boxSize {};
  long timeInc {}, time {};
  string positionFile {};
  bool fileOpen {false};
  ifstream reader;
  
  // For storing bead position and type data 
  vector<vector<double> >* position;
  vector<vector<int> >* boundaryCount;

  void initData();
  void deleteData();

public:
  
  // Destructor
  ~PositionReader();
  
  // Load a new position file
  bool open(string positionFile, int npoints, 
	    double lx, double ly, int timeInc);
  void close();

  // Return the position of a bead
  double getPosition(int index, int comp) const;
  double getUnwrappedPosition(int index, int comp) const;

  // Retrieve the bead position data in the next time frame
  // Return false if there is no next frame
  bool nextFrame();

  // Skip the bead position data in the next time frame
  bool skipFrame();
  
  // Get the time associated with the frame
  long getTime() const;

  // Check if there is a file open
  bool isOpen();
  
};

#endif
