// neighbour.h
// A code to find the nearest neighbours of every cell in the phase field model

#ifndef NEIGHBOUR_H
#define NEIGHBOUR_H

#ifdef __cplusplus
extern "C" {
#endif
  
#include "phase_field_model.h"
  
  typedef struct NeighbourAnalyser {
    int lx;
    int ly;
    double** totalField;
    int** indexField;
  } NeighbourAnalyser;
  
  typedef struct NeighbourList {
    int cellIndex; // Index of this cell
    int numOfNeighbours;
    int* neighbourIndex; // Store the list of neighbour cell indices
  } NeighbourList;
  
  NeighbourAnalyser* createNeighbourAnalyser(int lx, int ly);
  void deleteNeighbourAnalyser(NeighbourAnalyser* ana);
  void deleteNeighbourList(NeighbourList* list);
  NeighbourList** getNeighbourList(NeighbourAnalyser* ana,
				   PhaseFieldModel* model);
  int** getIndexField(NeighbourAnalyser* ana, PhaseFieldModel* model);
  
#ifdef __cplusplus
}
#endif
#endif
