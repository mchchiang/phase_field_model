// dump_overlap_field.c
// Dump the local overlap field for each cell

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"
#include "overlap.h"
#include "array.h"
#include "constant.h"

typedef struct OverlapFieldDump {
  Dump super; // Base struct must be the first element
  int cellIndex;
  bool overwrite;
  OverlapAnalyser* analyser;
} OverlapFieldDump;

void overlapFieldOutput(OverlapFieldDump* dump, PhaseFieldModel* model, 
			int step) {
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
  int clx = dump->analyser->cellLx;
  int cly = dump->analyser->cellLy;
  computeOverlap(dump->analyser, model, dump->cellIndex);
  for (int i = 0; i < clx; i++) {
    for (int j = 0; j < cly; j++) {
      fprintf(f, "%d %d %.5f\n", i, j, dump->analyser->overlapField[i][j]);
    }
    fprintf(f, "\n");
  }  
  fclose(f);
  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }
}

void deleteOverlapFieldDump(OverlapFieldDump* dump) {
  deleteOverlapAnalyser(dump->analyser);
  free(dump);
}

DumpFuncs overlapFieldDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &overlapFieldOutput,
   .destroy = (void (*)(Dump*)) &deleteOverlapFieldDump
  };

Dump* createOverlapFieldDump(char* filename, int clx, int cly, int cellIndex,
			     int printInc, bool overwrite) {
  OverlapFieldDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &overlapFieldDumpFuncs);
  dump->cellIndex = cellIndex;
  dump->overwrite = overwrite;
  dump->analyser = createOverlapAnalyser(clx, cly);
  return (Dump*) dump;
}
