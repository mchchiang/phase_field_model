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
  int nsteps = stoi(string(argv[1]), nullptr, 10);
  int ncells = stoi(string(argv[2]), nullptr, 10);
  double timeInc = stod(string(argv[3]), nullptr);
  PhaseFieldModel* model = new PhaseFieldModel(100, 100, 1);
  model->setDiffusionCoeff(0, 1.0);
  model->setIdealCellVolume(0, 120);
  model->setRegulateCoeff(0, 1.0);
  model->setVolumeCoeff(0, 1.0);
  model->setExclusionCoeff(0, 0, 1.0);
  model->setAdhesionCoeff(0, 0, 0.0);
  model->setMotility(0, 0.7);
  model->initCellLattice(ncells, 0, 30, 30);
  model->setDt(timeInc);
  model->run(nsteps);
  delete model;
}
