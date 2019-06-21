// neighbour.cpp

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <set>
using std::vector;
using std::set;

extern "C" {
#include "cell.h"
#include "array.h"
#include "phase_field_model.h"
#include "neighbour.h"

  NeighbourAnalyser* createNeighbourAnalyser(int lx, int ly) {
    NeighbourAnalyser* ana = (NeighbourAnalyser*) malloc(sizeof *ana);
    ana->lx = lx;
    ana->ly = ly;
    ana->totalField = create2DDoubleArray(lx, ly);
    ana->indexField = create2DIntArray(lx, ly);
    return ana;
  }
  
  void deleteNeighbourAnalyser(NeighbourAnalyser* ana) {
    free(ana->totalField);
    free(ana->indexField);
    free(ana);
  }
  
  void deleteNeighbourList(NeighbourList* list) {
    free(list->neighbourIndex);
    free(list);
  }
  
  NeighbourList** getNeighbourList(NeighbourAnalyser* ana,
				   PhaseFieldModel* model) {
    // Reset the fields
    for (int i = 0; i < ana->lx; i++) {
      for (int j = 0; j < ana->ly; j++) {
	ana->totalField[i][j] = 0.0;
	ana->indexField[i][j] = 0;
      }
    }
    
    // Map all individual cell fields to the total field
    int clx, cly, cx, cy, x, y, get;
    double phi, threshold;
    Cell* cell;
    int ncells = model->numOfCells;
    for (int i = 0; i < ncells; i++) {
      cell = model->cells[i];
      clx = cell->lx;
      cly = cell->ly;
      cx = cell->x;
      cy = cell->y;
      get = cell->getIndex;
      threshold = cell->incell;
      for (int j = 0; j < clx; j++) {
	for (int k = 0; k < cly; k++) {
	  x = iwrap(model, cx+j);
	  y = jwrap(model, cy+k);
	  phi = cell->field[get][j][k];
	  if (phi > ana->totalField[x][y] && phi > threshold) {
	    ana->totalField[x][y] = phi;
	    ana->indexField[x][y] = i+1; // Add one and use zero for substrate
	  }
	}
      }
    }
    
    // Find neighbours
    vector<set<int> > neighbours (ncells);
    int cellIndex, neighbourIndex;
    int nearest[2];
    for (int i = 0; i < ana->lx; i++) {
      for (int j = 0; j < ana->ly; j++) {
	cellIndex = ana->indexField[i][j];
	if (cellIndex == 0) continue; // Ignore substrate
	nearest[0] = ana->indexField[iup(model,i)][j];
	nearest[1] = ana->indexField[i][jup(model,j)];
	for (int k = 0; k < 2; k++) {
	  neighbourIndex = nearest[k];
	  if (cellIndex != neighbourIndex && neighbourIndex != 0) {
	    // Output cell index should start from zero
	    neighbours[cellIndex-1].insert(neighbourIndex-1);
	    neighbours[neighbourIndex-1].insert(cellIndex-1);
	  }
	}
      }
    }
    
    // Create the neighbour list
    NeighbourList** list = (NeighbourList**)
      malloc(sizeof *list * ncells);
    for (int i = 0; i < ncells; i++) {
      list[i] = (NeighbourList*) malloc(sizeof(NeighbourList));
      list[i]->cellIndex = i;
      list[i]->numOfNeighbours = neighbours[i].size();
      list[i]->neighbourIndex = create1DIntArray(list[i]->numOfNeighbours);
      int j = 0;
      for (set<int>::iterator it = neighbours[i].begin();
	   it != neighbours[i].end(); it++) {
	list[i]->neighbourIndex[j] = *it;
	j++;
      }
    }
    
    return list;
  }
}
