/*
 * CellGroup.cpp
 *
 *  Created on: 21 Dec, 2018
 *      Author: MichaelChiang
 */

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#include "CellGroup.hpp"
#include "Cell.hpp"

using std::vector;
using std::string;
using std::ostringstream;
using std::endl;

CellGroup::CellGroup(int _lx, int _ly) {
  lx = _lx;
  ly = _ly;
  vector<double> zeros = vector<double>(ly, 0.0);
  field = vector<vector<double> >(lx, zeros);
  cells = vector<Cell*>();
}

CellGroup::~CellGroup() {
  cells.clear();
}

void CellGroup::updateField() {
  // reset the field to zero
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      field[i][j] = 0.0;
    }
  }

  for (const Cell* cell : cells) {
    int cellLx = cell->getLx();
    int cellLy = cell->getLy();
    int cellX = cell->getX();
    int cellY = cell->getY();
    int x, y;
    for (int i = 0; i < cellLx; i++) {
      for (int j = 0; j < cellLy; j++) {
        x = iwrap(cellX + i);
        y = jwrap(cellY + j);
        field[x][y] += cell->getVolume(i, j);
      }
    }
  }
}

string CellGroup::printField() {
  ostringstream oss;
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      oss << field[i][j] << " ";
    }
    oss << endl;
  }
  return oss.str();
}

void CellGroup::addCell(Cell* cell) {
  cells.push_back(cell);
}

void CellGroup::removeCell(Cell* cell) {
  cells.erase(std::find(cells.begin(), cells.end(), cell));
}


int CellGroup::getNumOfCells() const {
  return static_cast<int>(cells.size());
}

int CellGroup::getLx() const {
  return lx;
}

int CellGroup::getLy() const {
  return ly;
}

double CellGroup::get(int i, int j) const {
  return field[i][j];
}

void CellGroup::set(int i, int j, double value) {
  field[i][j] = value;
}

int CellGroup::iwrap(int i) {
  int remainder = i % lx;
  if (remainder >= 0) {
    return remainder;
  }
  return lx + remainder;
}

int CellGroup::jwrap(int j) {
  int remainder = j % ly;
  if (remainder >= 0) {
    return remainder;
  }
  return ly + remainder;
}
