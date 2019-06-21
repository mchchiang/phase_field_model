// dump_bulk_cm.c
// Dump the bulk cm of all the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"
#include "array.h"

typedef struct BulkCMDump {
  Dump super; // Base struct must be the first element
} BulkCMDump;

void bulkCMOutput(BulkCMDump* dump, PhaseFieldModel* model, int step) {
  FILE* f;
  f = fopen(dump->super.filename, "a");
  // Compute the bulk cm
  Cell* cell;
  int clx, cly, get;
  int lx = model->lx;
  int ly = model->ly;
  double mass = 0.0;
  double totalMass = 0.0;
  double xcm = 0.0;
  double ycm = 0.0;
  for (int i = 0; i < model->numOfCells; i++) {
    cell = model->cells[i];
    clx = cell->lx;
    cly = cell->ly;
    get = cell->getIndex;
    mass = 0.0;
    for (int j = 0; j < clx; j++) {
      for (int k = 0; k < cly; k++) {
	mass += cell->field[get][j][k];
      }
    }
    totalMass += mass;
    xcm += ((model->cellXCM[i]+model->cellXBoundCount[i]*lx)*mass);
    ycm += ((model->cellYCM[i]+model->cellYBoundCount[i]*ly)*mass);
  }
  xcm /= totalMass;
  ycm /= totalMass;
  
  fprintf(f, "%d %.5f %.5f\n", step, xcm, ycm);
  fclose(f);
}

void deleteBulkCMDump(BulkCMDump* dump) {
  free(dump);
}

DumpFuncs bulkCMDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &bulkCMOutput,
   .destroy = (void (*)(Dump*)) &deleteBulkCMDump
  };

Dump* createBulkCMDump(char* filename, int printInc) {
  BulkCMDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &bulkCMDumpFuncs);
  return (Dump*) dump;
}
