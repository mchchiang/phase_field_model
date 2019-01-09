/*
 * PhaseFieldModel.cpp
 *
 *  Created on: 21 Dec, 2018
 *      Author: MichaelChiang
 */

#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
//#include <omp.h>

#include "Cell.hpp"
#include "CellGroup.hpp"
#include "PhaseFieldModel.hpp"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::ofstream;

PhaseFieldModel::PhaseFieldModel(int _lx, int _ly, int numOfCellGroups) {
  lx = _lx;
  ly = _ly;

  int numOfCellGroupsSq = numOfCellGroups * numOfCellGroups;
  cellGroups = vector<CellGroup*>(numOfCellGroups, new CellGroup(lx, ly));
  D = vector<double>(numOfCellGroups, 0.0);
  alpha = vector<double>(numOfCellGroups, 0.0);
  idealCellVolume = vector<double>(numOfCellGroups, 0.0);
  gamma = vector<double>(numOfCellGroups, 0.0);
  beta = vector<double>(numOfCellGroupsSq, 0.0);
  eta = vector<double>(numOfCellGroupsSq, 0.0);
  cellLx = vector<int>(numOfCellGroups, 0);
  cellLy = vector<int>(numOfCellGroups, 0);
  motility = vector<double>(numOfCellGroups, 0.0);
}

PhaseFieldModel::~PhaseFieldModel() {
  for (int i = 0; i < cellGroups.size(); i++) {
    delete cellGroups[i];
  }
  cellGroups.clear();

  for (int i = 0; i < cells.size(); i++) {
    delete cells[i];
  }
  cells.clear();
}

void PhaseFieldModel::initCellLattice(int numOfCells, int type,
                                      int cx, int cy) {
  double sqrtNumOfCells {sqrt(numOfCells)};
  int dx {static_cast<int>(floor(lx/sqrtNumOfCells))};
  int dy {static_cast<int>(floor(ly/sqrtNumOfCells))};
  int nx {lx/dx};
  int ny {ly/dy};
  cellLx[type] = cx;
  cellLy[type] = cy;

  int labx, laby, cellx, celly;
  int x0 {(cx-dx)/2};
  int y0 {(cy-dy)/2};

  for (int i = 0; i < numOfCells; i++) {
    labx = dx*(i/ny);
    laby = dy*(i%ny);
    cellx = labx-x0;
    celly = laby-y0;
    CellGroup* group = cellGroups[type];
    Cell* cell = new Cell(cellx, celly, cellLx[type], cellLy[type], type);
    cell->initOnes(x0, y0, dx, dy);
    cell->setTheta(0.0);
    group->addCell(cell);
    cells.push_back(cell);
  }
}

void PhaseFieldModel::run(int nsteps) {

  for (int i = 0; i < nsteps; i++) {
    for (int j = 0; j < cellGroups.size(); j++) {
      updateCellGroupVolume(cellGroups[j]);
    }

#pragma omp parallel for schedule(auto)
    for (int j = 0; j < cells.size(); j++) {
      updateCellVolume(cells[j]);
      updateCellField(cells[j]);
      updateCellCM(cells[j]);
      updateVelocity(cells[j]);
    }

    output(i);
  }
}

void PhaseFieldModel::output(int step) {
  if (step % 1000 == 0){
    cout << "Step " << step << endl;
    // Output lattice
    ofstream writer ("output.dat");
    for (int i = 0; i < lx; i++) {
      for (int j = 0; j < ly; j++) {
        writer << cellGroups[0]->get(i, j) << " ";
      }
      writer << endl;
    }
    writer.close();
  }
}

void PhaseFieldModel::updateCellGroupVolume(CellGroup* group) {
  group->updateField();
}

void PhaseFieldModel::updateCellVolume(Cell* cell) {
  cell->updateTotalVolume();
}

void PhaseFieldModel::updateCellCM(Cell* cell) {
  cell->updateCM();
}

void PhaseFieldModel::updateVelocity(Cell* cell) {
  cell->updateVelocity();
}

void PhaseFieldModel::updateCellField(Cell* cell) {
  // Apply fixed (Dirichlet) boundary conditions (u = 0 at boundaries)
  // i and j are coordinates in the cell's reference frame
  cell->startUpdateCellField();
  int type = cell->getCellType();
  int cx = cellLx[type];
  int cy = cellLy[type];
  for (int i = 1; i < cx - 1; i++) {
    for (int j = 1; j < cy - 1; j++) {
      cell->set(i, j, cell->get(i, j) + dt * (
          singleCellInteractions(cell, i, j) +
          cellCellInteractions(cell, i, j) +
          cellSubstrateInteractions(cell, i, j)));
    }
  }
  cell->endUpdateCellField();
}

double PhaseFieldModel::singleCellInteractions(Cell* cell, int i, int j) {
  double u = cell->get(i, j);
  double px = cell->getPx();
  double py = cell->getPy();
  int type = cell->getCellType();
  return D[type] * centralDiff(i, j, cell) - motility[type] * (
      px * forwardDiff(i, j, 0, cell) - py * forwardDiff(i, j, 1, cell)) +
      u * (1 - u) * (u - 0.5 + getVolumeCoeff(type)
          * (getIdealCellVolume(type) - cell->getTotalVolume()));
}

double PhaseFieldModel::cellCellInteractions(Cell* cell, int i, int j) {
  double u = cell->get(i, j);
  int type = cell->getCellType();
  int cellX = cell->getX();
  int cellY = cell->getY();
  int x = iwrap(cellX + i);
  int y = jwrap(cellY + j);

  CellGroup* group;

  // Compute excluded volume effect
  double exclusion = getExclusionCoeff(type, type) * cell->getVolume(i, j);
  for (int l = 0; l < cellGroups.size(); l++) {
    group = cellGroups[l];
    exclusion -= getExclusionCoeff(l, type) * group->get(x, y);
  }

  // Compute adhesion effect
  double volumeFieldCentralDiff = centralDiff(i, j, cell->getVolumeField());
  double adhesion = -getAdhesionCoeff(type, type) * volumeFieldCentralDiff;
  for (int l = 0; l < cellGroups.size(); l++) {
    group = cellGroups[l];
    adhesion += getAdhesionCoeff(l, type) * centralDiff(x, y, group);
  }

  // Compute regularisation effect
  double regularisation = getRegulateCoeff(type) * volumeFieldCentralDiff;

  return u * (1 - u) * (exclusion + adhesion + regularisation);

}

double PhaseFieldModel::cellSubstrateInteractions(Cell* cell, int i, int j) {
  return 0.0; // Not considered at the moment
}

// Accessor methods
int PhaseFieldModel::getLx() const {
  return lx;
}

int PhaseFieldModel::getLy() const {
  return ly;
}

void PhaseFieldModel::setDt(double _dt) {
  if (_dt > 0.0) {
    this->dt = _dt;
  }
}

double PhaseFieldModel::getDt() const {
  return dt;
}

int PhaseFieldModel::getCellLx(int type) const {
  return cellLx[type];
}

int PhaseFieldModel::getCellLy(int type) const {
  return cellLy[type];
}

int PhaseFieldModel::getNumOfCells() const {
  return static_cast<int>(cells.size());
}

int PhaseFieldModel::getNumOfCellGroups() const {
  return static_cast<int>(cellGroups.size());
}

CellGroup* PhaseFieldModel::getCellGroup(int type) {
  return cellGroups[type];
}

void PhaseFieldModel::setIdealCellVolume(int type, double value) {
  if (value >= 0.0) {
    idealCellVolume[type] = value;
  }
}

double PhaseFieldModel::getIdealCellVolume(int type) const {
  return idealCellVolume[type];
}

void PhaseFieldModel::setVolumeCoeff(int type, double value) {
  if (value >= 0.0) {
    alpha[type] = value;
  }
}

double PhaseFieldModel::getVolumeCoeff(int type) const {
  return alpha[type];
}

void PhaseFieldModel::setExclusionCoeff(int type1, int type2, double value) {
  if (value >= 0.0) {
    beta[getTypeIndex(type1, type2)] = value;
    beta[getTypeIndex(type2, type1)] = value;
  }
}

double PhaseFieldModel::getExclusionCoeff(int type1, int type2) const {
  return beta[getTypeIndex(type1, type2)];
}

void PhaseFieldModel::setAdhesionCoeff(int type1, int type2, double value) {
  if (value >= 0.0) {
    eta[getTypeIndex(type1, type2)] = value;
    eta[getTypeIndex(type2, type1)] = value;
  }
}

double PhaseFieldModel::getAdhesionCoeff(int type1, int type2) const {
  return eta[getTypeIndex(type1, type2)];
}

void PhaseFieldModel::setRegulateCoeff(int type, double value) {
  if (value >= 0.0) {
    gamma[type] = value;
  }
}

double PhaseFieldModel::getRegulateCoeff(int type) const {
  return gamma[type];
}

void PhaseFieldModel::setDiffusionCoeff(int type, double value) {
  if (value >= 0.0) {
    D[type] = value;
  }
}

double PhaseFieldModel::getDiffusionCoeff(int type) const {
  return D[type];
}

void PhaseFieldModel::setMotility(int type, double value) {
  if (value >= 0.0) {
    motility[type] = value;
  }
}

double PhaseFieldModel::getMotility(int type) const {
  return motility[type];
}

// Implement periodic boundary conditions
int PhaseFieldModel::iwrap(int i) {
  int remainder = i % lx;
  if (remainder >= 0) {
    return remainder;
  }
  return lx + remainder;
}

int PhaseFieldModel::jwrap(int j) {
  int remainder = j % ly;
  if (remainder >= 0) {
    return remainder;
  }
  return ly + remainder;
}

int PhaseFieldModel::iup(int i) {
  if (i + 1 >= lx) {
    return 0;
  }
  return i + 1;
}

int PhaseFieldModel::idown(int i) {
  if (i - 1 < 0) {
    return lx - 1;
  }
  return i - 1;
}

int PhaseFieldModel::jup(int j) {
  if (j + 1 >= ly) {
    return 0;
  }
  return j + 1;
}

int PhaseFieldModel::jdown(int j) {
  if (j - 1 < 0) {
    return ly - 1;
  }
  return j - 1;
}

double PhaseFieldModel::forwardDiff(int i, int j, int comp, Field2D* field){
  if (comp == 0) {
    return 0.5*(field->get(iup(i), j) - field->get(idown(i), j));
  }
  return 0.5*(field->get(i, jup(j)) - field->get(i, jdown(j)));
}

double PhaseFieldModel::centralDiff(int i, int j, Field2D* field) {
  return field->get(iup(i), j) + field->get(idown(i), j) + field->get(i, jup(j))
      + field->get(i, jdown(j)) - 4.0 * (field->get(i, j));
}

int PhaseFieldModel::getTypeIndex(int type1, int type2) const {
  return type1*getNumOfCellGroups()+type2;
}
