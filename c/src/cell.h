// cell.h

#ifndef CELL_H
#define CELL_H

#include "mtwister.h"

#define INCELL 0.4
#define CMSHIFT 2.0

typedef struct {
  double** field[2];
  double** volumeField;
  int setIndex;
  int getIndex;
  int lx; // x size of lattice
  int ly; // y size of lattice
  int x; // x pos of the cell relative to main lattice
  int y; // y pos of the cell relative to main lattice
  int type;
  double vx; // x velocity of the cell
  double vy; // y velocity of the cell
  double theta;
  double rotateDiff;
  MTRand random; 
  double xcm; // x centre of mass in cell's frame
  double ycm; // y centre of mass in cell's frame
  double deltaXCM;
  double deltaYCM;
  double volume; // total volume of the cell
} Cell;

void initCell(Cell* cell, int x, int y, int lx, int ly, int type, int seed);
void deleteCell(Cell* cell);
void initFieldSquareOnes(Cell* cell, int x0, int y0, int dx, int dy);
void updateCM(Cell* cell);
void updateVolume(Cell* cell);
void updateVelocity(Cell* cell);
void shiftCoordinates(Cell* cell, int xShift, int yShift);
void calculateCM(Cell* cell, double* xcm, double* ycm);
void startUpdateCellField(Cell* cell);
void endUpdateCellField(Cell* cell);

#endif
