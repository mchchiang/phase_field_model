/*
 * Cell.cpp
 *
 *  Created on: 21 Dec, 2018
 *      Author: MichaelChiang
 */

#include <cmath>
#include <vector>
#include <iostream>

#include "Cell.hpp"

using std::vector;
using std::cout;
using std::endl;

Cell::Cell (int _x, int _y, int _lx, int _ly, int type) {
  x = _x;
  y = _y;
  lx = _lx;
  ly = _ly;
  xcm = lx / 2.0;
  ycm = ly / 2.0;
  deltaXCM = 0.0;
  deltaYCM = 0.0;
  cellType = type;

  // Create the cell field to zero
  vector<double> zeros {vector<double>(ly, 0.0)};
  for (int k = 0; k < 2; k++) {
    cellField[k] = vector<vector<double> >(lx, zeros);
  }

  this->setField = 1;
  this->getField = 1;
  this->volumeField = VolumeField(this);
}

void Cell::initSquareCell(int dx, int dy) {
  int xStart = static_cast<int>(lx/2.0 - dx/2.0);
  int xEnd = xStart + dx;
  int yStart = static_cast<int>(ly/2.0 - dy/2.0);
  int yEnd = yStart + dy;
  for (int k = 0; k < 2; k++){
    for (int i = xStart; i < xEnd; i++) {
      for (int j = yStart; j < yEnd; j++) {
        cellField[k][i][j] = 1.0;
      }
    }
  }
  vector<double> cm = calculateCM();
  xcm = cm[0];
  ycm = cm[1];
}

void Cell::initCell(const vector<vector<double> >& matrix) {
  if (static_cast<int>(matrix.size()) == lx &&
      static_cast<int>(matrix[0].size()) == ly) {
    for (int i = 0; i < lx; i++) {
      for (int j = 0; j < ly; j++) {
        cellField[0][i][j] = matrix[i][j];
        cellField[1][i][j] = matrix[i][j];
      }
    }
  }
  vector<double> cm = calculateCM();
  xcm = cm[0];
  ycm = cm[1];
}

void Cell::updateTotalVolume() {
  volume = calculateTotalVolume();
}

double Cell::calculateTotalVolume() {
  double sum = 0.0;
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      sum += volumeField.get(i, j);
    }
  }
  return sum;
}

vector<double> Cell::calculateCM() {
  double xavg = 0.0;
  double yavg = 0.0;
  int count = 0;
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      if (cellField[getField][i][j] > INCELL) {
        xavg += i;
        yavg += j;
        count++;
      }
    }
  }

  if (count > 0) {
    xavg /= static_cast<double>(count);
    yavg /= static_cast<double>(count);
  } else {
    xavg = 0.0;
    yavg = 0.0;
  }
  return {xavg, yavg};
}

void Cell::updateCM() {
  double oldXCM = xcm;
  double oldYCM = ycm;

  vector<double> cm = calculateCM();
  xcm = cm[0];
  ycm = cm[1];

  deltaXCM += (xcm - oldXCM);
  deltaYCM += (ycm - oldYCM);

  if (fabs(deltaXCM) > CMSHIFT || fabs(deltaYCM) > CMSHIFT) {
    int xShift = static_cast<int>(floor(deltaXCM));
    int yShift = static_cast<int>(floor(deltaYCM));
    shiftCoordinates(xShift, yShift);
    deltaXCM -= xShift;
    deltaYCM -= yShift;
    cm = calculateCM();
    xcm = cm[0];
    ycm = cm[1];
  }
}

void Cell::shiftCoordinates(int xShift, int yShift) {
  startUpdateCellField();
  int xStart, xEnd, yStart, yEnd;
  int zeroXStart, zeroXEnd, zeroYStart, zeroYEnd;
  if (xShift >= 0) {
    xStart = xShift;
    xEnd = lx;
    zeroXStart = lx - xShift;
    zeroXEnd = lx;
  } else {
    xStart = 0;
    xEnd = lx + xShift;
    zeroXStart = 0;
    zeroXEnd = -xShift;
  }
  if (yShift >= 0) {
    yStart = yShift;
    yEnd = ly;
    zeroYStart = ly - yShift;
    zeroYEnd = ly;
  } else {
    yStart = 0;
    yEnd = ly + yShift;
    zeroYStart = 0;
    zeroYEnd = -yShift;
  }

  // Set empty cells to zero
  for (int i = zeroXStart; i < zeroXEnd; i++) {
    for (int j = 0; j < ly; j++) {
      set(i, j, 0.0);
    }
  }
  for (int i = zeroXEnd; i < lx; i++) {
    for (int j = zeroYStart; j < zeroYEnd; j++) {
      set(i, j, 0.0);
    }
  }

  // Shift cells
  for (int i = xStart; i < xEnd; i++) {
    for (int j = yStart; j < yEnd; j++) {
      set(i - xShift, j - yShift, get(i, j));
    }
  }

  // Update end coordinate
  x += xShift;
  y += yShift;
  endUpdateCellField();
}

double Cell::getXCM() const {
  return xcm;
}

double Cell::getYCM() const {
  return ycm;
}

int Cell::getX() const {
  return x;
}

int Cell::getY() const {
  return y;
}

double Cell::getVolume(int i, int j) const {
  return volumeField.get(i, j);
}

double Cell::getTotalVolume() const {
  return volume;
}

Field2D* Cell::getVolumeField() {
  return &volumeField;
}

int Cell::getCellType() const {
  return cellType;
}

int Cell::getLx() const {
  return lx;
}

int Cell::getLy() const {
  return ly;
}

double Cell::get(int i, int j) const {
  return cellField[getField][i][j];
}

void Cell::set(int i, int j, double value) {
  cellField[setField][i][j] = value;
}

void Cell::startUpdateCellField() {
  setField = (getField == 1 ? 0 : 1);
}

void Cell::endUpdateCellField() {
  getField = setField;
}
