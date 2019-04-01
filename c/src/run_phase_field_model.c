// run_phase_field.c

#include <stdio.h>
#include <stdlib.h>
#include "phase_field_model.h"

int main (int argc, char* argv[]) {
  int argi = 0;

  if (argc != 7) {
    printf("Usage: RunPhaseFieldModel [L] [ncells] [cellVol] [cellL] " \
	   "[nsteps] [timeInc]\n");
    return 1;
  }

  int L = atoi(argv[++argi]);
  int ncells = atoi(argv[++argi]);
  double cellVol = atof(argv[++argi]);
  int cellL = atoi(argv[++argi]);
  int nsteps = atoi(argv[++argi]);
  double timeInc = atof(argv[++argi]);

  printf("Initialising model ...\n");
  
  PhaseFieldModel* model = (PhaseFieldModel*) malloc(sizeof(PhaseFieldModel));
  initModel(model, L, L, 1, ncells);
  model->diffusion[0] = 0.4;
  model->idealCellVolume[0] = cellVol;
  model->gamma[0] = 4.0;
  model->alpha[0] = 0.0025;
  model->beta[getTypeIndex(model, 0, 0)] = 1.0;
  model->eta[getTypeIndex(model, 0, 0)] = 3.0;
  model->motility[0] = 0.0;
  model->cellLx[0] = cellL;
  model->cellLy[0] = cellL;
  int x = 0, y = 0;
  for (int i = 0; i < ncells; i++) {
    initSquareCell(model, i, x, y, 10, 10, 0);
    y += 8;
    if (y >= L) {
      y = x / 8 + 4;
      x += 8;
    }
  }

  printf("Done initialisation.\n");
  
  model->dt = timeInc;
  run(model, 30000);
  model->motility[0] = 0.5;
  run(model, nsteps);
  deleteModel(model);
}
