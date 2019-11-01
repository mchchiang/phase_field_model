// cell.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cell.h"
#include "array.h"
#include "random.h"
#include "constant.h"

Cell* createCell(int x, int y, int lx, int ly,
		 double dr, double incell, unsigned long seed) {
  // Allocate memory for a cell
  Cell* cell = malloc(sizeof *cell);
  
  cell->x = x;
  cell->y = y;
  cell->lx = lx;
  cell->ly = ly;
  cell->incell = incell;

  // Allocate memory and initialise the field to zero
  for (int i = 0; i < 2; i++) {
    cell->field[i] = create2DDoubleArray(cell->lx, cell->ly);
  }
  cell->setIndex = 1;
  cell->getIndex = 1;
  cell->rotateDiff = dr;
  cell->xcm = cell->lx/2.0;
  cell->ycm = cell->ly/2.0;
  cell->deltaXCM = 0.0;
  cell->deltaYCM = 0.0;
  cell->volume = 0.0;
  cell->random = createRandom(seed);
  // pick a random 2D direction
  cell->theta = randDouble(cell->random)*2.0*PF_PI;
  cell->vx = cos(cell->theta);
  cell->vy = sin(cell->theta);
  cell->v = 1.0;
  cell->gyration = create1DDoubleArray(3);
  return cell;
}

void deleteCell(Cell* cell) {
  for (int i = 0; i < 2; i++) {
    free(cell->field[i]);
  }
  deleteRandom(cell->random);
  free(cell->gyration);
  free(cell);
}

void initField(Cell* cell, double** lattice) {
  for (int i = 0; i < cell->lx; i++) {
    for (int j = 0; j < cell->ly; j++) {
      cell->field[0][i][j] = lattice[i][j];
      cell->field[1][i][j] = lattice[i][j];
    }
  }
  calculateCM(cell, &cell->xcm, &cell->ycm);
  updateVolume(cell);
}

void initFieldSquare(Cell* cell, int x0, int y0, int dx, int dy, double phi0) {
  for (int i = 0; i < dx && x0+i < cell->lx-1; i++) {
    for (int j = 0; j < dy && y0+j < cell->ly-1; j++) {
      cell->field[0][i+x0][j+y0] = phi0;
      cell->field[1][i+x0][j+y0] = phi0;
    }
  }
  calculateCM(cell, &cell->xcm, &cell->ycm);
  updateVolume(cell);
}

void calculateCM(Cell* cell, double* xcm, double* ycm) {
  double xavg = 0.0;
  double yavg = 0.0;
  double mass = 0.0;
  double phi;
  int get = cell->getIndex;
  for (int i = 0; i < cell->lx; i++) {
    for (int j = 0; j < cell->ly; j++) {
      phi = cell->field[get][i][j];
      xavg += (phi * (i+0.5)); // Use the centre of a lattice element
      yavg += (phi * (j+0.5));
      mass += phi;
    }
  }
  if (mass > 0.0) {
    xavg /= mass;
    yavg /= mass;
  } else {
    xavg = 0.0;
    yavg = 0.0;
  }
  *xcm = xavg;
  *ycm = yavg;
}

void updateCM(Cell* cell) {
  double oldXCM = cell->xcm;
  double oldYCM = cell->ycm;
  calculateCM(cell, &cell->xcm, &cell->ycm);
  cell->deltaXCM += (cell->xcm - oldXCM);
  cell->deltaYCM += (cell->ycm - oldYCM);

  if (fabs(cell->deltaXCM) > PF_CMSHIFT || fabs(cell->deltaYCM) > PF_CMSHIFT) {
    int xShift = (int)(round(cell->deltaXCM));
    int yShift = (int)(round(cell->deltaYCM));
    shiftCoordinates(cell, xShift, yShift);
    cell->deltaXCM -= xShift;
    cell->deltaYCM -= yShift;
    cell->xcm -= xShift;
    cell->ycm -= yShift;
  }
}

void updateVolume(Cell* cell) {
  double totalVolume = 0.0;
  int get = cell->getIndex;
  double phi;
  for (int i = 0; i < cell->lx; i++) {
    for (int j = 0; j < cell->ly; j++) {
      phi = cell->field[get][i][j];
      totalVolume += phi*phi;
    }
  }
  cell->volume = totalVolume;
}

void updateVelocity(Cell* cell, double dt) {
  cell->theta += sqrt(6.0 * cell->rotateDiff * dt) *
    (randDouble(cell->random)*2.0-1.0);
  cell->vx = cos(cell->theta);
  cell->vy = sin(cell->theta);
}

void shiftCoordinates(Cell* cell, int xShift, int yShift) {
  startUpdateCellField(cell);
  int set = cell->setIndex;
  int get = cell->getIndex;
  int xStart, xEnd, yStart, yEnd; // these are in the old frame
  int zeroXStart, zeroXEnd, zeroYStart, zeroYEnd; // these are in the new frame
  if (xShift >= 0) {
    xStart = xShift;
    xEnd = cell->lx;
    zeroXStart = cell->lx - xShift;
    zeroXEnd = cell->lx;
  } else {
    xStart = 0;
    xEnd = cell->lx + xShift;
    zeroXStart = 0;
    zeroXEnd = -xShift;
  }
  if (yShift >= 0) {
    yStart = yShift;
    yEnd = cell->ly;
    zeroYStart = cell->ly - yShift;
    zeroYEnd = cell->ly;
  } else {
    yStart = 0;
    yEnd = cell->ly + yShift;
    zeroYStart = 0;
    zeroYEnd = -yShift;
  }

  // Set empty cells to zero
  for (int i = zeroXStart; i < zeroXEnd; i++) {
    for (int j = 0; j < cell->ly; j++) {
      cell->field[set][i][j] = 0.0;
    }
  }

  if (zeroYStart < zeroYEnd) {
    for (int i = 0; i < cell->lx; i++) {
      for (int j = zeroYStart; j < zeroYEnd; j++) {
	cell->field[set][i][j] = 0.0;
      }
    }
  }

  // Shift cells
  for (int i = xStart; i < xEnd; i++) {
    for (int j = yStart; j < yEnd; j++) {
      cell->field[set][i-xShift][j-yShift] = cell->field[get][i][j];
    }
  }

  // Update end coordinate
  cell->x += xShift;
  cell->y += yShift;
  endUpdateCellField(cell);
}

void startUpdateCellField(Cell* cell) {
  cell->setIndex = (cell->getIndex == 1 ? 0 : 1);
}

void endUpdateCellField(Cell* cell) {
  cell->getIndex = cell->setIndex;
}
