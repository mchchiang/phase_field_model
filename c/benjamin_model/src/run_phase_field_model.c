// run_phase_field.c

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "dump.h"
#include "phase_field_model.h"

int main (int argc, char* argv[]) {
  
  if (argc != 2) {
    printf("Usage: RunPhaseFieldModel param_file [dump_file]\n");
    return 1;
  }

  int argi = 0;
  char* filename = argv[++argi];
  
  printf("Reading parameters ...\n");

  FILE* file = fopen(filename, "r");
  if (file == NULL) {
    printf("Error in opening parameter file!\n");
    return 1;
  }

  char line [80], dumpMode [80];
  char cmFile [DIR_SIZE], shapeFile [DIR_SIZE], dumpFile [DIR_SIZE];
  double phi0, M, R, kappa, alpha, mu, Dr, epsilon, dt, v;
  int cellLx, cellLy, lx, ly, nequil, nsteps, ncells, seed;
  int nparams = 0;
  int ndumps = 0;
  int nedumps = 0;
  int maxDumps = 50;
  int printInc, override, cellIndex;
  Dump** equilDumps = malloc(sizeof *equilDumps * maxDumps);
  Dump** dumps = malloc(sizeof *dumps * maxDumps);

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
    
    // Read dumps
    if (sscanf(line, "dump_cm %d %d %s %s", 
	       &printInc, &override, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = createCMDump(dumpFile, printInc, override);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createCMDump(dumpFile, printInc, override);
	ndumps++; 
      }
    }
    if (sscanf(line, "dump_gyr %d %d %s %s",
	       &printInc, &override, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = 
	  createGyrationDump(dumpFile, printInc, override);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createGyrationDump(dumpFile, printInc, override);
	ndumps++;
      }
    }
    if (sscanf(line, "dump_field %d %d %s %s",
	       &printInc, &override, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = createFieldDump(dumpFile, printInc, override);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createFieldDump(dumpFile, printInc, override);
	ndumps++;
      }
    }
    if (sscanf(line, "dump_cell_field %d %d %d %s %s",
	       &cellIndex, &printInc, &override, dumpMode, dumpFile) == 5) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] =
	  createCellFieldDump(dumpFile, cellIndex, printInc, override);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] =
	  createCellFieldDump(dumpFile, cellIndex, printInc, override);
	ndumps++;
      }
    }
  }
  equilDumps = realloc(equilDumps, sizeof *equilDumps * nedumps);
  dumps = realloc(dumps, sizeof *dumps * ndumps);

  fclose(file);
  
  if (nparams != 20) {
    printf("Not enough parameters specified!\n");
    return 1;
  }  
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
  printf("shape_file = %s\n", shapeFile);  

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
  model->ndumps = nedumps;
  model->dumps = equilDumps;
  
  initCellsFromFile(model, cmFile, shapeFile, seed);
  
  printf("Done initialisation.\n");
  
  model->dt = dt;

  printf("Doing equilibration run ...\n");
  run(model, nequil);

  model->motility = v;
  model->ndumps = ndumps;
  model->dumps = dumps;

  printf("Doing main simulation run ...\n");
  run(model, nsteps);

  deleteModel(model);

  for (int i = 0; i < nedumps; i++) {
    deleteDump(equilDumps[i]);
  }
  free(equilDumps);

  for (int i = 0; i < ndumps; i++) {
    deleteDump(dumps[i]);
  }
  free(dumps);
}
