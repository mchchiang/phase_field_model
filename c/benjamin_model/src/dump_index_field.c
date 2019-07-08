// dump_index_field.c
// Dump the index field

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"
#include "array.h"
#include "neighbour.h"

typedef struct IndexFieldDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
} IndexFieldDump;

void indexFieldOutput(IndexFieldDump* dump, PhaseFieldModel* model, int step) {
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
  int** field = getIndexField(NULL, model);
  for (int i = 0; i < model->lx; i++) {
    for (int j = 0; j < model->ly; j++) {
      fprintf(f, "%d %d %d\n", i, j, field[i][j]);
    }
    fprintf(f, "\n");
  }
  fclose(f);
  free(field);
  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }
}

void deleteIndexFieldDump(IndexFieldDump* dump) {
  free(dump);
}

DumpFuncs indexFieldDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &indexFieldOutput,
   .destroy = (void (*)(Dump*)) &deleteIndexFieldDump
  };

Dump* createIndexFieldDump(char* filename, int printInc, bool overwrite) {
  IndexFieldDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &indexFieldDumpFuncs);
  dump->overwrite = overwrite;
  return (Dump*) dump;
}
