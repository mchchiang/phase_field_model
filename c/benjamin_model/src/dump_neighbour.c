// dump_neighbours.c
// Dump the neighbour list of all the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"
#include "neighbour.h"

typedef struct NeighbourDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
  NeighbourAnalyser* analyser;
} NeighbourDump;

void neighbourOutput(NeighbourDump* dump, PhaseFieldModel* model, int step) {
  char tmpfile [PF_DIR_SIZE];
  FILE* f;
  if (dump->overwrite) {
    strcpy(tmpfile, dump->super.filename);
    strcat(tmpfile, ".tmp");
    f = fopen(tmpfile, "w");
  } else {
    f = fopen(dump->super.filename, "a");
  }

  // Find neighbours
  fprintf(f, "Cells: %d\n", model->numOfCells);
  fprintf(f, "Timestep: %d\n", step);
  NeighbourList** list = getNeighbourList(dump->analyser, model);
  for (int i = 0; i < model->numOfCells; i++) {
    for (int j = 0; j < list[i]->numOfNeighbours; j++) {
      fprintf(f, "%d ", list[i]->neighbourIndex[j]);
    }
    fprintf(f,"\n");
  }
  fclose(f);

  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }
  for (int i = 0; i < model->numOfCells; i++) {
    deleteNeighbourList(list[i]);
  }
  free(list);
}

void deleteNeighbourDump(NeighbourDump* dump) {
  deleteNeighbourAnalyser(dump->analyser);
  free(dump);
}

DumpFuncs neighbourDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &neighbourOutput,
   .destroy = (void (*)(Dump*)) &deleteNeighbourDump
  };

Dump* createNeighbourDump(char* filename, int lx, int ly,
			  int printInc, bool overwrite) {
  NeighbourDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &neighbourDumpFuncs);
  dump->overwrite = overwrite;
  dump->analyser = createNeighbourAnalyser(lx, ly);
  return (Dump*) dump;
}
