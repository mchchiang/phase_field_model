/*
 * RunPhaseFieldModel.cpp
 *
 *  Created on: 9 Jan, 2019
 *      Author: MichaelChiang
 */

#include <string>
#include "PhaseFieldModel.hpp"

using std::string;

int main (int argc, char* argv[]) {
  int argi = 0;
  int l = stoi(string(argv[++argi]), nullptr, 10);
  int ncells = stoi(string(argv[++argi]), nullptr, 10);
  double cellVol = stod(string(argv[++argi]), nullptr);
  int cellL = stoi(string(argv[++argi]), nullptr, 10);
  int nsteps = stoi(string(argv[++argi]), nullptr, 10);
  double timeInc = stod(string(argv[++argi]), nullptr);
  PhaseFieldModel* model = new PhaseFieldModel(l, l, 1);
  model->setDiffusionCoeff(0, 1.0);
  model->setIdealCellVolume(0, cellVol);
  model->setRegulateCoeff(0, 10.0);
  model->setVolumeCoeff(0, 1.0);
  model->setExclusionCoeff(0, 0, 10.0);
  model->setAdhesionCoeff(0, 0, 0.0);
  model->setMotility(0, 1.0);
  model->initCellLattice(ncells, 0, cellL, cellL);
  model->setDt(timeInc);
  model->run(nsteps);
  delete model;
}
