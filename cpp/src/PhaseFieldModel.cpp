/*
 * PhaseFieldModel.cpp
 *
 *  Created on: 21 Dec, 2018
 *      Author: MichaelChiang
 */

#include <cmath>
#include <vector>
#include <random>
#include <iostream>

#include "Cell.hpp"
#include "CellGroup.hpp"
#include "PhaseFieldModel.hpp"

using std::vector;
using std::cout;
using std::endl;

PhaseFieldModel::PhaseFieldModel(int lx, int ly, int numOfCellGroups) {
  this->lx = lx;
  this->ly = ly;

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
}

PhaseFieldModel::~PhaseFieldModel() {
  CellGroup* group;
  while (cellGroups.size() > 0) {
    group = cellGroups[0];
    delete group;
  }

  Cell* cell;
  while (cells.size() > 0) {
    cell = cells[0];
    delete cell;
  }
}

void PhaseFieldModel::initSquareCellLattice(int x0, int y0, int xlen, int ylen,
                                            int cx, int cy, int numOfCells,
                                            int type){
  double volumePerCell = xlen*ylen/static_cast<double>(numOfCells);
  int dx = static_cast<int>(floor(sqrt(volumePerCell)));
  int dy = dx;
  int x = x0;
  int y = y0;

  for (int i = 0; i < numOfCells; i++){
    /*std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> randInt(0,1);
    type = randInt(mt);*/
    CellGroup* group = cellGroups[type];
    cellLx[type] = cx;
    cellLy[type] = cy;
    double idealVolume = idealCellVolume[type];
    int cellLen = static_cast<int>(floor(sqrt(idealVolume)));
    Cell* cell = new Cell(x, y, cellLx[type], cellLy[type], type);
    cell->initSquareCell(cellLen, cellLen);
    group->addCell(cell);
    cells.push_back(cell);
    x += dx;
    if (x > x0 + xlen){
      x = x0;
      y += dy;
    }
  }
}

void PhaseFieldModel::run(int nsteps) {

  for (size_t i {}; i < nsteps; i++) {

    for (size_t j {}; j < cellGroups.size(); j++){
      updateCellGroupVolume(cellGroups[j]);
    }

    for (size_t j {}; j < cells.size(); j++){
      updateCellVolume(cells[j]);
      updateCellField(cells[j]);
      updateCellCM(cells[j]);
    }

    output(i);
  }
}

void PhaseFieldModel::output(int step){
  if (step % 1000 == 0){
    cout << "Step " << step << endl;
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

void PhaseFieldModel::updateCellField(Cell* cell) {
  // Apply fixed (Dirichlet) boundary conditions (u = 0 at boundaries)
  // i and j are coordinates in the cell's reference frame
  cell->startUpdateCellField();
  int type = cell->getCellType();
  int cx = cellLx[type];
  int cy = cellLy[type];
  for (size_t i = 1; i < cx - 1; i++) {
    for (size_t j = 1; j < cy - 1; j++) {
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
  int type = cell->getCellType();
  return D[type] * centralDiff(i, j, *cell) +
      0.05*(forwardDiff(i, j, 0, *cell) + forwardDiff(i, j, 1, *cell)) +
      u * (1 - u) * (u - 0.5 + getVolumeCoeff(type)
          * (getIdealCellVolume(type) - cell->getTotalVolume()));
}

double PhaseFieldModel::cellCellInteractions(Cell *cell, int i, int j) {
  double u = cell->get(i, j);
  int type = cell->getCellType();
  int cellX = cell->getX();
  int cellY = cell->getY();
  int x = iwrap(cellX + i);
  int y = jwrap(cellY + j);

  CellGroup* group;

  // Compute excluded volume effect
  double exclusion = getExclusionCoeff(type, type) * cell->getVolume(i, j);
  for (size_t l = 0; l < cellGroups.size(); l++) {
    group = cellGroups[l];
    exclusion -= getExclusionCoeff(l, type) * group->get(x, y);
  }

  // Compute adhesion effect
  double volumeFieldCentralDiff = centralDiff(i, j, *(cell->getVolumeField()));
  double adhesion = -getAdhesionCoeff(type, type) * volumeFieldCentralDiff;
  for (size_t l = 0; l < cellGroups.size(); l++) {
    group = cellGroups[l];
    adhesion += getAdhesionCoeff(l, type) * centralDiff(x, y, *group);
  }

  // Compute regularisation effect
  double regularisation = getRegulateCoeff(type) * volumeFieldCentralDiff;

  return u * (1 - u) * (exclusion + adhesion + regularisation);
}

double PhaseFieldModel::cellSubstrateInteractions(Cell* cell, int i, int j) {
  return 0.0; // Not considered at the moment
}

// Accessor methods
int PhaseFieldModel::getLx() {
  return lx;
}

int PhaseFieldModel::getLy() {
  return ly;
}

void PhaseFieldModel::setDt(double dt){
  if (dt > 0.0){
    this->dt = dt;
  }
}

double PhaseFieldModel::getDt(){
  return dt;
}

int PhaseFieldModel::getCellLx(int type){
  return cellLx[type];
}

int PhaseFieldModel::getCellLy(int type){
  return cellLy[type];
}

int PhaseFieldModel::getNumOfCells() {
  return static_cast<int>(cells.size());
}

int PhaseFieldModel::getNumOfCellGroups() {
  return static_cast<int>(cellGroups.size());
}

CellGroup* PhaseFieldModel::getCellGroup(int type){
  return cellGroups[type];
}

void PhaseFieldModel::setIdealCellVolume(int type, double value){
  if (value >= 0.0){
    idealCellVolume[type] = value;
  }
}

double PhaseFieldModel::getIdealCellVolume(int type) {
  return idealCellVolume[type];
}

void PhaseFieldModel::setVolumeCoeff(int type, double value){
  if (value >= 0.0){
    alpha[type] = value;
  }
}

double PhaseFieldModel::getVolumeCoeff(int type) {
  return alpha[type];
}

void PhaseFieldModel::setExclusionCoeff(int type1, int type2, double value){
  if (value >= 0.0){
    beta[getTypeIndex(type1, type2)] = value;
    beta[getTypeIndex(type2, type1)] = value;
  }
}

double PhaseFieldModel::getExclusionCoeff(int type1, int type2) {
  return beta[getTypeIndex(type1, type2)];
}

void PhaseFieldModel::setAdhesionCoeff(int type1, int type2, double value){
  if (value >= 0.0){
    eta[getTypeIndex(type1, type2)] = value;
    eta[getTypeIndex(type2, type1)] = value;
  }
}

double PhaseFieldModel::getAdhesionCoeff(int type1, int type2) {
  return eta[getTypeIndex(type1, type2)];
}

void PhaseFieldModel::setRegulateCoeff(int type, double value) {
  if (value >= 0.0){
    gamma[type] = value;
  }
}

double PhaseFieldModel::getRegulateCoeff(int type) {
  return gamma[type];
}

void PhaseFieldModel::setDiffusionCoeff(int type, double value){
  if (value >= 0.0){
    D[type] = value;
  }
}

double PhaseFieldModel::getDiffusionCoeff(int type) {
  return D[type];
}

// Implement periodic boundary conditions
int PhaseFieldModel::iwrap(int i) {
  int remainder = i % lx;
  if (remainder >= 0)
    return remainder;
  return lx + remainder;
}

int PhaseFieldModel::jwrap(int j) {
  int remainder = j % ly;
  if (remainder >= 0)
    return remainder;
  return ly + remainder;
}

int PhaseFieldModel::iup(int i) {
  if (i + 1 >= lx)
    return 0;
  return i + 1;
}

int PhaseFieldModel::idown(int i) {
  if (i - 1 < 0)
    return lx - 1;
  return i - 1;
}

int PhaseFieldModel::jup(int j) {
  if (j + 1 >= ly)
    return 0;
  return j + 1;
}

int PhaseFieldModel::jdown(int j) {
  if (j - 1 < 0)
    return ly - 1;
  return j - 1;
}

double PhaseFieldModel::forwardDiff(int i, int j, int comp, Field2D& field){
  if (comp == 0){
    return 0.5*(field.get(iup(i), j) - field.get(idown(i), j));
  }
  return 0.5*(field.get(i, jup(j)) - field.get(i, jdown(j)));
}

double PhaseFieldModel::centralDiff(int i, int j, Field2D& field) {
  return field.get(iup(i), j) + field.get(idown(i), j) + field.get(i, jup(j))
      + field.get(i, jdown(j)) - 4.0 * (field.get(i, j));
}

int PhaseFieldModel::getTypeIndex(int type1, int type2){
  return type1*getNumOfCellGroups()+type2;
}

int main (int argc, char* argv[]){
  PhaseFieldModel* model = new PhaseFieldModel(80, 80, 1);
  model->setDiffusionCoeff(0, 1.0);
  model->setIdealCellVolume(0, 25);
  model->setRegulateCoeff(0, 1.0);
  model->setVolumeCoeff(0, 1.0);
  model->setExclusionCoeff(0, 0, 1.0);
  model->setAdhesionCoeff(0, 0, 0.3);
  model->initSquareCellLattice(10, 10, 70, 70, 20, 20, 80, 0);
  model->run(10000);
}
