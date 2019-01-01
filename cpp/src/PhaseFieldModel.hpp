/*
 * PhaseFieldModel.hpp
 *
 *  Created on: 21 Dec, 2018
 *      Author: MichaelChiang
 */

#ifndef PHASEFIELDMODEL_HPP_
#define PHASEFIELDMODEL_HPP_

#include <vector>

#include "Cell.hpp"
#include "CellGroup.hpp"

class PhaseFieldModel {

private:
  int lx {}, ly {};
  double dt {0.01};

  std::vector<double> idealCellVolume; // Ideal cell volume
  std::vector<double> alpha; // Growth coefficient
  std::vector<double> beta; // Excluded volume coefficient
  std::vector<double> gamma; // Regularisation for adhesion
  std::vector<double> eta; // Adhesion coefficient
  std::vector<double> D; // Surface-tension-like coefficient
  std::vector<int> cellLx;
  std::vector<int> cellLy;

  std::vector<CellGroup*> cellGroups;
  std::vector<Cell*> cells;

protected:
  int iwrap(int i);
  int jwrap(int j);
  int iup(int i);
  int idown(int i);
  int jup(int j);
  int jdown(int j);
  int getTypeIndex(int type1, int type2);
  double centralDiff(int i, int j, Field2D* field);
  double forwardDiff(int i, int j, int comp, Field2D* field);

public:
  PhaseFieldModel(int lx, int ly, int numOfCellGroups);
  ~PhaseFieldModel();

  void initSquareCellLattice(int x0, int y0, int xlen, int ylen,
                            int cx, int cy, int numOfCells, int type);

  void run(int nsteps);
  void updateCellGroupVolume(CellGroup* group);
  void updateCellVolume(Cell* cell);
  void updateCellField(Cell* cell);
  void updateCellCM(Cell* cell);
  void output(int step);

  double singleCellInteractions(Cell* cell, int i, int j);
  double cellCellInteractions(Cell* cell, int i, int j);
  double cellSubstrateInteractions(Cell* cell, int i, int j);

  // Accessor methods
  int getLx();
  int getLy();
  void setDt(double dt);
  double getDt();
  int getCellLx(int type);
  int getCellLy(int type);
  int getNumOfCells();
  int getNumOfCellGroups();
  CellGroup* getCellGroup(int type);
  void setIdealCellVolume(int type, double value);
  double getIdealCellVolume(int type);
  void setVolumeCoeff(int type, double value);
  double getVolumeCoeff(int type);
  void setExclusionCoeff(int type1, int type2, double value);
  double getExclusionCoeff(int type1, int type2);
  void setAdhesionCoeff(int type1, int type2, double value);
  double getAdhesionCoeff(int type1, int type2);
  void setRegulateCoeff(int type, double value);
  double getRegulateCoeff(int type);
  void setDiffusionCoeff(int type, double value);
  double getDiffusionCoeff(int type);

};

#endif /* PHASEFIELDMODEL_HPP_ */
