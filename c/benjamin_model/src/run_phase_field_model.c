// run_phase_field.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "phase_field_model.h"

int main (int argc, char* argv[]) {
  
  if (argc != 2) {
    printf("Usage: RunPhaseFieldModel [param_file]\n");
    return 1;
  }

  int argi = 0;
  char* filename = argv[++argi];
  
  printf("Reading parameters ...\n");

  FILE* file = fopen(filename, "r");
  if (file == NULL) {
    printf("Error in opening file!\n");
    return 1;
  }

  char line [80], cmFile[500], shapeFile[500];
  double phi0, M, R, kappa, alpha, mu, Dr, epsilon, dt, v;
  int cellLx, cellLy, lx, ly, nequil, nsteps, ncells, seed;
  int nparams = 0;
  while (fgets(line, sizeof(line), file) != NULL) {
    nparams += sscanf(line, "alpha = %lf", &alpha);
    nparams += sscanf(line, "phi0 = %lf", &phi0);
    nparams += sscanf(line, "M = %lf", &M);
    nparams += sscanf(line, "R = %lf", &R);
    nparams += sscanf(line, "kappa = %lf", &kappa);
    nparams += sscanf(line, "mu = %lf", &mu);
    nparams += sscanf(line, "epsilon = %lf", &epsilon);
    nparams += sscanf(line, "cellLx = %d", &cellLx);
    nparams += sscanf(line, "cellLy = %d", &cellLy);
    nparams += sscanf(line, "lx = %d", &lx);
    nparams += sscanf(line, "ly = %d", &ly);
    nparams += sscanf(line, "nsteps = %d", &nsteps);
    nparams += sscanf(line, "nequil = %d", &nequil);
    nparams += sscanf(line, "ncells = %d", &ncells);
    nparams += sscanf(line, "dt = %lf", &dt);
    nparams += sscanf(line, "Dr = %lf", &Dr);
    nparams += sscanf(line, "v = %lf", &v);
    nparams += sscanf(line, "cm_file = %s", cmFile);
    nparams += sscanf(line, "shape_file = %s", shapeFile);
    nparams += sscanf(line, "seed = %d", &seed);
  }

  if (nparams != 20) {
    printf("Not enough parameters specified!\n");
    return 1;
  } else {
    printf("Read parameters:\n");
    printf("lx = %d\n", lx);
    printf("ly = %d\n", ly);
    printf("cellLx = %d\n", cellLx);
    printf("cellLy = %d\n", cellLy);
    printf("ncells = %d\n", ncells);
    printf("phi0 = %.5f\n", phi0);
    printf("mu = %.5f\n", mu);
    printf("alpha = %.5f\n", alpha);
    printf("kappa = %.5f\n", kappa);
    printf("epsilon = %.5f\n", epsilon);
    printf("Dr = %.5f\n", Dr);
    printf("dt = %.5f\n", dt);
    printf("v = %.5f\n", v);
    printf("M = %.5f\n", M);
    printf("R = %.5f\n", R);
    printf("seed = %d\n", seed);
    printf("cm_file = %s\n", cmFile);
    printf("shape)file = %s\n", shapeFile);
  }

  printf("Initialising model ...\n");
  
  PhaseFieldModel* model = createModel(lx, ly, ncells);
  model->phi0 = phi0;
  model->M = M;
  model->piR2phi02 = M_PI * R * R * phi0 * phi0;
  model->kappa = kappa;
  model->alpha = alpha;
  model->mu = mu;
  model->Dr = Dr;
  model->epsilon = epsilon;  
  model->motility = 0.0;
  model->cellLx = cellLx;
  model->cellLy = cellLy;

  initCellsFromFile(model, cmFile, shapeFile, seed);
  // Init cell in the middle
  // initSquareCell(model, 0, 0, 0, (int)R, (int)R);
  /*  int x = 0, y = 0;
  for (int i = 0; i < ncells; i++) {
    initSquareCell(model, i, x, y, (int)R, (int)R);
    y += (int) (R*1.4);
    if (y >= L) {
      y = (int) (x / R + (R*0.5));
      x += (int) (R*1.2);
    }
    }*/

  printf("Done initialisation.\n");
  
  model->dt = dt;
  run(model, nequil);
  model->mu = mu;
  model->motility = v;
  run(model, nsteps);
  deleteModel(model);
}
