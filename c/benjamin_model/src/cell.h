// cell.h

#ifndef CELL_H
#define CELL_H

#include "random.h"

typedef struct {
  double** field[2]; // phase field
  int setIndex;
  int getIndex;
  int lx; // x size of lattice
  int ly; // y size of lattice
  int x; // x pos of the cell relative to main lattice
  int y; // y pos of the cell relative to main lattice
  int type;
  double v; // Speed of the cell
  double vx; // x velocity of the cell
  double vy; // y velocity of the cell
  double theta;
  double rotateDiff;
  Random* random; 
  double xcm; // x centre of mass in cell's frame
  double ycm; // y centre of mass in cell's frame
  double deltaXCM;
  double deltaYCM;
  double volume; // total volume of the cell
  double incell;
  double* gyration; // gyration tensor of the cell
} Cell;

Cell* createCell(int x, int y, int lx, int ly, double dr,
		 double incell, unsigned long seed);
void deleteCell(Cell* cell);
void initFieldSquare(Cell* cell, int x0, int y0, int dx, int dy, double phi0);
void initField(Cell* cell, double** field);
void updateCM(Cell* cell);
void updateVolume(Cell* cell);
void updateVelocity(Cell* cell, double dt);
void shiftCoordinates(Cell* cell, int xShift, int yShift);
void calculateCM(Cell* cell, double* xcm, double* ycm);
void startUpdateCellField(Cell* cell);
void endUpdateCellField(Cell* cell);

#endif
