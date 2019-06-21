// dump_field.c
// Dump the entire field

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"
#include "array.h"

typedef struct FieldDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
} FieldDump;

void fieldOutput(FieldDump* dump, PhaseFieldModel* model, int step) {
  char tmpfile [PF_DIR_SIZE];
  if (dump->overwrite) {
    strcpy(tmpfile, dump->super.filename);
    strcat(tmpfile, ".tmp");
  } else {
    strcpy(tmpfile, dump->super.filename);
    char suffix [80];
    sprintf(suffix, ".%d", step);
    strcat(tmpfile, suffix);
  }
  FILE* f = fopen(tmpfile, "w");
  // Sum over all the phase fields to find the total field
  Cell* cell;
  int get, clx, cly, cx, cy, x, y;
  double phi;
  double** field = create2DDoubleArray(model->lx, model->ly);
  for (int i = 0; i < model->numOfCells; i++) {
    cell = model->cells[i];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    get = cell->getIndex;
    for (int j = 0; j < clx; j++) {
      for (int k = 0; k < cly; k++) {
	x = iwrap(model, cx+j);
	y = jwrap(model, cy+k);
	phi = cell->field[get][j][k];
	field[x][y] += phi;
      }
    }
  }
  for (int i = 0; i < model->lx; i++) {
    for (int j = 0; j < model->ly; j++) {
      fprintf(f, "%d %d %.5f\n", i, j, field[i][j]);
    }
    fprintf(f, "\n");
  }
  fclose(f);
  free(field);
  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }
}

void deleteFieldDump(FieldDump* dump) {
  free(dump);
}

DumpFuncs fieldDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &fieldOutput,
   .destroy = (void (*)(Dump*)) &deleteFieldDump
  };

Dump* createFieldDump(char* filename, int printInc, bool overwrite) {
  FieldDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &fieldDumpFuncs);
  dump->overwrite = overwrite;
  return (Dump*) dump;
}
