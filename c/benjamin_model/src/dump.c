// dump.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"

void setDump(Dump* dump, void* derived, char* filename, int printInc,
	     DumpFuncs* funcs) {
  dump->derived = derived;
  dump->filename = filename;
  dump->printInc = printInc;
  dump->funcs = funcs;
}

void deleteDump(Dump* dump) {
  dump->funcs->delete(dump);
  free(dump);
}

void dumpOutput(Dump* dump, PhaseFieldModel* model, int step) {
  if (step % dump->printInc == 0) {
    dump->funcs->output(dump, model, step);
  }
}

// ************************************************************************* //

// Implementations of specific dumpers

// Dump the entire field

typedef struct FieldDump {
  Dump super; // Base struct must be the first element
} FieldDump;

void fieldOutput(FieldDump* dump, PhaseFieldModel* model, int step) {
  char tmpfile [80];
  strcpy(tmpfile, dump->super.filename);
  strcat(tmpfile, ".tmp");
  FILE* f = fopen(tmpfile, "w");
  for (int i = 0; i < model->lx; i++) {
    for (int j = 0; j < model->ly; j++) {
      fprintf(f, "%d %d %.5f\n", i, j, sqrt(model->totalField[i][j]));
    }
    fprintf(f, "\n");
  }  
  fclose(f);
  rename(tmpfile, dump->super.filename);
}

void deleteFieldDump(FieldDump* dump) {}

DumpFuncs fieldDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &fieldOutput,
   .delete = (void (*)(Dump*)) &deleteFieldDump
  };

Dump* createFieldDump(char* filename, int printInc) {
  FieldDump* dump = (FieldDump*) malloc(sizeof(FieldDump));
  setDump(&dump->super, dump, filename, printInc, &fieldDumpFuncs);
  return (Dump*) dump;
}

// ************************************************************************* //

// Dump the field for a specific cell

typedef struct CellFieldDump {
  Dump super;
  int cellIndex;
} CellFieldDump;

void cellFieldOutput(CellFieldDump* dump, PhaseFieldModel* model, int step) {
  char tmpfile [80];
  strcpy(tmpfile, dump->super.filename);
  strcat(tmpfile, ".tmp");
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
  rename(tmpfile, dump->super.filename);
}

void deleteCellFieldDump(CellFieldDump* dump) {}

DumpFuncs cellFieldDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &cellFieldOutput,
   .delete = (void (*)(Dump*)) &deleteCellFieldDump
  };

Dump* createCellFieldDump(char* filename, int cellIndex, int printInc) {
  CellFieldDump* dump = (CellFieldDump*) malloc(sizeof(CellFieldDump));
  setDump(&dump->super, dump, filename, printInc, &cellFieldDumpFuncs);
  dump->cellIndex = cellIndex;
  return (Dump*) dump;
}

// ************************************************************************* //

// Dump the cm of the cells

typedef struct CMDump {
  Dump super; // Base struct must be the first element
} CMDump;

void cmOutput(CMDump* dump, PhaseFieldModel* model, int step) {
  char tmpfile [80];
  strcpy(tmpfile, dump->super.filename);
  strcat(tmpfile, ".tmp");
  FILE* f = fopen(tmpfile, "w");
  Cell* cell;
  double x, y, cx, cy;
  for (int i = 0; i < model->numOfCells; i++) {
    cell = model->cells[i];
    cx = cell->x;
    cy = cell->y;
    x = iwrap(model, cx+cell->xcm);
    y = jwrap(model, cy+cell->ycm);
    fprintf(f, "%d %.5f %.5f\n", i, x, y);
  }
  fclose(f);
  rename(tmpfile, dump->super.filename);
}

void deleteCMDump(CMDump* dump) {}

DumpFuncs cmDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &cmOutput,
   .delete = (void (*)(Dump*)) &deleteCMDump
  };

Dump* createCMDump(char* filename, int printInc) {
  CMDump* dump = (CMDump*) malloc(sizeof(CMDump));
  setDump(&dump->super, dump, filename, printInc, &cmDumpFuncs);
  return (Dump*) dump;
}
