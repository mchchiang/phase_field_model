/*
 * RunPhaseFieldModel.cpp
 *
 *  Created on: 9 Jan, 2019
 *      Author: MichaelChiang
 */

#include <string>
#include <iostream>
#include "PhaseFieldModel.hpp"

using std::string;
using std::cout;
using std::endl;

int main (int argc, char* argv[]) {

  int argi {};

  if (argc != 7) {
    cout << "Usage: RunPhaseFieldModel [L] [ncells] [cellVol] [cellL] "
         << "[nsteps] [timeInc]" << endl;
    return 1;
  }

  int L {stoi(string(argv[++argi]), nullptr, 10)};
  int ncells {stoi(string(argv[++argi]), nullptr, 10)};
  double cellVol {stod(string(argv[++argi]), nullptr)};
  int cellL {stoi(string(argv[++argi]), nullptr, 10)};
  int nsteps {stoi(string(argv[++argi]), nullptr, 10)};
  double timeInc {stod(string(argv[++argi]), nullptr)};

  PhaseFieldModel* model = new PhaseFieldModel(L, L, 1);
  model->setDiffusionCoeff(0, 0.4);
  model->setIdealCellVolume(0, cellVol);
  model->setRegulateCoeff(0, 4.0);
  model->setVolumeCoeff(0, 0.0025);
  model->setExclusionCoeff(0, 0, 1.0);
  model->setAdhesionCoeff(0, 0, 3.0);
  model->setMotility(0, 0.0);

  model->initSquareCell(40, 44, 10, 10, cellL, cellL, 0);
  model->initSquareCell(40, 36, 10, 10, cellL, cellL, 0);
  model->initSquareCell(30, 24, 10, 10, cellL, cellL, 0);
  model->initSquareCell(30, 26, 10, 10, cellL, cellL, 0);
  model->initSquareCell(20, 74, 10, 10, cellL, cellL, 0);
  model->initSquareCell(20, 76, 10, 10, cellL, cellL, 0);
  model->initSquareCell(10, 64, 10, 10, cellL, cellL, 0);
  model->initSquareCell(10, 66, 10, 10, cellL, cellL, 0);
  model->initSquareCell(50, 14, 10, 10, cellL, cellL, 0);
  model->initSquareCell(50, 16, 10, 10, cellL, cellL, 0);
  model->initSquareCell(45, 44, 10, 10, cellL, cellL, 0);
  model->initSquareCell(45, 36, 10, 10, cellL, cellL, 0);
  model->initSquareCell(35, 24, 10, 10, cellL, cellL, 0);
  model->initSquareCell(35, 26, 10, 10, cellL, cellL, 0);
  model->initSquareCell(25, 74, 10, 10, cellL, cellL, 0);
  model->initSquareCell(25, 76, 10, 10, cellL, cellL, 0);
  model->initSquareCell(15, 64, 10, 10, cellL, cellL, 0);
  model->initSquareCell(15, 66, 10, 10, cellL, cellL, 0);
  model->initSquareCell(65, 14, 10, 10, cellL, cellL, 0);
  model->initSquareCell(65, 16, 10, 10, cellL, cellL, 0);
  model->initSquareCell(75, 44, 10, 10, cellL, cellL, 0);
  model->initSquareCell(75, 36, 10, 10, cellL, cellL, 0);
  model->initSquareCell(85, 24, 10, 10, cellL, cellL, 0);
  model->initSquareCell(85, 26, 10, 10, cellL, cellL, 0);
  model->initSquareCell(95, 74, 10, 10, cellL, cellL, 0);
  model->initSquareCell(95, 76, 10, 10, cellL, cellL, 0);
  model->initSquareCell(5, 64, 10, 10, cellL, cellL, 0);
  model->initSquareCell(5, 66, 10, 10, cellL, cellL, 0);
  model->initSquareCell(70, 14, 10, 10, cellL, cellL, 0);
  model->initSquareCell(70, 16, 10, 10, cellL, cellL, 0);

  model->setDt(timeInc);
  model->run(20000);
  model->setMotility(0, 1.0);
  model->run(nsteps);
  delete model;
}
