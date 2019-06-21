// dump_cell_field.c
// Dump the field for a specific cell

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"

typedef struct CellFieldDump {
  Dump super;
  int cellIndex;
  bool overwrite;
} CellFieldDump;

void cellFieldOutput(CellFieldDump* dump, PhaseFieldModel* model, int step) {
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
  Cell* cell = model->cells[dump->cellIndex];
  int get = cell->getIndex;
  for (int i = 0; i < model->cellLx; i++) {
    for (int j = 0; j < model->cellLy; j++) {
      fprintf(f, "%d %d %.5f\n", i, j, cell->field[get][i][j]);
    }
    fprintf(f, "\n");
  }  
  fclose(f);
  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }
}

void deleteCellFieldDump(CellFieldDump* dump) {
  free(dump);
}

DumpFuncs cellFieldDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &cellFieldOutput,
   .destroy = (void (*)(Dump*)) &deleteCellFieldDump
  };

Dump* createCellFieldDump(char* filename, int cellIndex,
			  int printInc, bool overwrite) {
  CellFieldDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &cellFieldDumpFuncs);
  dump->cellIndex = cellIndex;
  dump->overwrite = overwrite;
  return (Dump*) dump;
}
