// run_phase_field.c

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
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
  double phi0 = -1.0;
  double M, R, kappa, alpha, mu, Dr, epsilon, dt, v;
  int cellLx = -1;
  int cellLy = -1;
  int lx, ly, nequil, nsteps, ncells;
  unsigned long seed;
  int nparams = 0;
  int ndumps = 0;
  int nedumps = 0;
  int maxDumps = 50;
  int printInc, overwrite, cellIndex;
  int fieldScale, kernelLength, sgolayDegree, sgolayLength;
  double kernelSigma;
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
    nparams += sscanf(line, "seed = %ld", &seed);
    
    // Read dumps
    // CM dump
    if (sscanf(line, "dump_cm %d %d %s %s", 
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = createCMDump(dumpFile, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createCMDump(dumpFile, printInc, overwrite);
	ndumps++; 
      }
    }
    // Bulk CM dump
    if (sscanf(line, "dump_bulk_cm %d %s %s",
	       &printInc, dumpMode, dumpFile) == 3) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = createBulkCMDump(dumpFile, printInc);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createBulkCMDump(dumpFile, printInc);
	ndumps++; 
      }
    }
    // Gyration dump
    if (sscanf(line, "dump_gyr %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = 
	  createGyrationDump(dumpFile, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createGyrationDump(dumpFile, printInc, overwrite);
	ndumps++;
      }
    }
    // Total field dump
    if (sscanf(line, "dump_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = createFieldDump(dumpFile, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createFieldDump(dumpFile, printInc, overwrite);
	ndumps++;
      }
    }
    // Individual cell field dump
    if (sscanf(line, "dump_cell_field %d %d %d %s %s",
	       &cellIndex, &printInc, &overwrite, dumpMode, dumpFile) == 5) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] =
	  createCellFieldDump(dumpFile, cellIndex, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] =
	  createCellFieldDump(dumpFile, cellIndex, printInc, overwrite);
	ndumps++;
      }
    }
    // Shape dump
    if (sscanf(line, "dump_shape %d %d %lf %d %d %d %d %s %s",
	       &fieldScale, &kernelLength, &kernelSigma,
	       &sgolayDegree, &sgolayLength,
	       &printInc, &overwrite, dumpMode, dumpFile) == 9) {
      // Only created the dump when cell field size and phi0 are known,
      // as they are needed for creating the shape analysers
      if (cellLx > 0 && cellLy > 0 && phi0 >= 0.0) {
	if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	  equilDumps[nedumps] =
	    createShapeDump(dumpFile, fieldScale, cellLx, cellLy,
			    kernelLength, kernelSigma, sgolayDegree,
			    sgolayLength, phi0/2.0, printInc, overwrite);
	  nedumps++;
	} else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	  dumps[ndumps] =
	    createShapeDump(dumpFile, fieldScale, cellLx, cellLy,
			    kernelLength, kernelSigma, sgolayDegree,
			    sgolayLength, phi0/2.0, printInc, overwrite);
	  ndumps++;
	}
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
  printf("seed = %ld\n", seed);
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

  double start, end, duration;
  printf("Doing equilibration run ...\n");
  start = omp_get_wtime();
  run(model, nequil);
  end = omp_get_wtime();
  duration = end-start;
  printf("Time taken (sec): %.5f\n", duration);
  printf("\n");

  model->motility = v;
  model->ndumps = ndumps;
  model->dumps = dumps;

  printf("Doing main simulation run ...\n");
  start = omp_get_wtime();
  run(model, nsteps);
  end = omp_get_wtime();
  duration = end-start;
  printf("Time taken (sec): %.5f\n", duration);
  
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
