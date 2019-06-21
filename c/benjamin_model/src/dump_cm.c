// dump_cm.c
// Dump the cm of the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"

typedef struct CMDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
} CMDump;

void cmOutput(CMDump* dump, PhaseFieldModel* model, int step) {
  char tmpfile [PF_DIR_SIZE];
  FILE* f;
  if (dump->overwrite) {
    strcpy(tmpfile, dump->super.filename);
    strcat(tmpfile, ".tmp");
    f = fopen(tmpfile, "w");
  } else {
    f = fopen(dump->super.filename, "a");
  }
  Cell* cell;
  int ix, iy;
  double x, y, cx, cy, wx, wy;
  fprintf(f, "Cells: %d\n", model->numOfCells);
  fprintf(f, "Timestep: %d\n", step);
  for (int i = 0; i < model->numOfCells; i++) {
    cell = model->cells[i];
    cx = cell->x;
    cy = cell->y;
    x = cx+cell->xcm;
    y = cy+cell->ycm;
    ix = (int) floor(x / model->lx);
    iy = (int) floor(y / model->ly);
    wx = x - ix * model->lx;
    wy = y - iy * model->ly;
    fprintf(f, "%.5f %.5f %d %d\n", wx, wy, ix, iy);
  }
  fclose(f);
  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }
}

void deleteCMDump(CMDump* dump) {
  free(dump);
}

DumpFuncs cmDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &cmOutput,
   .destroy = (void (*)(Dump*)) &deleteCMDump
  };

Dump* createCMDump(char* filename, int printInc, bool overwrite) {
  CMDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &cmDumpFuncs);
  dump->overwrite = overwrite;
  return (Dump*) dump;
}
