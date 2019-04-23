// dump_gyration.c
// Dump the gyration tensor components of the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"

typedef struct GyrationDump {
  Dump super; // Base struct must be the first element
  bool override;
} GyrationDump;

void gyrationOutput(GyrationDump* dump, PhaseFieldModel* model, int step) {
  char tmpfile [DUMP_FILE_SIZE];
  FILE* f;
  if (dump->override) {
    strcpy(tmpfile, dump->super.filename);
    strcat(tmpfile, ".tmp");
    f = fopen(tmpfile, "w");
  } else {
    f = fopen(dump->super.filename, "a");
  }
  Cell* cell;
  double gxx, gyy, gxy;
  fprintf(f, "Cells: %d\n", model->numOfCells);
  fprintf(f, "Timestep: %d\n", step);
  for (int i = 0; i < model->numOfCells; i++) {
    cell = model->cells[i];
    gxx = cell->gyration[0];
    gyy = cell->gyration[1];
    gxy = cell->gyration[2];
    fprintf(f, "%.5f %.5f %.5f\n", gxx, gyy, gxy);
  }
  fclose(f);
  if (dump->override) {
    rename(tmpfile, dump->super.filename);
  }
}

void deleteGyrationDump(GyrationDump* dump) {}

DumpFuncs gyrationDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &gyrationOutput,
   .delete = (void (*)(Dump*)) &deleteGyrationDump
  };

Dump* createGyrationDump(char* filename, int printInc, bool override) {
  GyrationDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &gyrationDumpFuncs);
  dump->override = override;
  return (Dump*) dump;
}
