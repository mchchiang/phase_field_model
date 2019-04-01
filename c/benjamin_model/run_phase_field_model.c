// run_phase_field.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "phase_field_model.h"

int main (int argc, char* argv[]) {
  int argi = 0;

  if (argc != 6) {
    printf("Usage: RunPhaseFieldModel [L] [ncells] [cellL] " \
	   "[nsteps] [timeInc]\n");
    return 1;
  }

  int L = atoi(argv[++argi]);
  int ncells = atoi(argv[++argi]);
  int cellL = atoi(argv[++argi]);
  int nsteps = atoi(argv[++argi]);
  double timeInc = atof(argv[++argi]);

  printf("Initialising model ...\n");
  
  PhaseFieldModel* model = (PhaseFieldModel*) malloc(sizeof(PhaseFieldModel));
  initModel(model, L, L, ncells);
  int R = 12;
  model->phi0 = 2.0;
  model->M = 0.1;
  model->piR2 = M_PI * R * R;
  model->kappa = 1.6;
  model->alpha = model->kappa/2.0;
  model->mu = 1.0;
  model->Dr = 0.0001;
  model->epsilon = 0.3;  
  model->motility = 0.0;
  model->cellLx = cellL;
  model->cellLy = cellL;
  int x = 0, y = 0;
  for (int i = 0; i < ncells; i++) {
    initSquareCell(model, i, x, y, R, R);
    y += (int) (R*1.4);
    if (y >= L) {
      y = (int) (x / R + (R*0.5));
      x += (int) (R*1.2);
    }
  }

  printf("Done initialisation.\n");
  
  model->dt = timeInc;
  run(model, 30000);
  model->mu = 10.0;
  model->motility = 0.5;
  run(model, nsteps);
  deleteModel(model);
}
