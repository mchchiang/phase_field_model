// dump_overlap.c
// Dump the local average overlap for each cell

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

typedef struct OverlapDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
  int nthreads; // Use multiple threads to process the images
  OverlapAnalyser** analysers; // One analyser per thread
} OverlapDump;

void overlapOutput(OverlapDump* dump, PhaseFieldModel* model, int step) {
  char tmpfile [PF_DIR_SIZE];
  FILE* f;
  if (dump->overwrite) {
    strcpy(tmpfile, dump->super.filename);
    strcat(tmpfile, ".tmp");
    f = fopen(tmpfile, "w");
  } else {
    f = fopen(dump->super.filename, "a");
  }
  // Compute the average local overlap for each cell
  int j;
  int id = 0;
  int ncells = model->numOfCells;
  OverlapAnalyser* analyser;
  double* overlapAvg = create1DDoubleArray(ncells);
#pragma omp parallel default(none) \
  shared(ncells, dump, model, overlapAvg) \
  private (j, id, analyser)
  {
#ifdef _OPENMP
    id = omp_get_thread_num();
#endif
    analyser = dump->analysers[id];
#pragma omp for schedule(static)
    for (j = 0; j < ncells; j++) {
      overlapAvg[j] = computeOverlap(analyser, model, j);
    }
  }
  // Output the average overlap data to file
  fprintf(f, "Cells: %d\n", ncells);
  fprintf(f, "Timestep: %d\n", step);
  for (int i = 0; i < ncells; i++) {
    fprintf(f, "%.5e\n", overlapAvg[i]);
  }
  fclose(f);
  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }

  // Clean up resources
  free(overlapAvg);
}

void deleteOverlapDump(OverlapDump* dump) {
  for (int i = 0; i < dump->nthreads; i++) {
    deleteOverlapAnalyser(dump->analysers[i]);
  }
  free(dump->analysers);
  free(dump);
}

DumpFuncs overlapDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &overlapOutput,
   .destroy = (void (*)(Dump*)) &deleteOverlapDump
  };

Dump* createOverlapDump(char* filename, int clx, int cly,
			int printInc, bool overwrite) {
  OverlapDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &overlapDumpFuncs);
  dump->overwrite = overwrite;
  // Initialise the overlap analysers
#ifdef _OPENMP
  dump->nthreads = omp_get_max_threads();
#else
  dump->nthreads = 1;
#endif
  dump->analysers = malloc(sizeof *dump->analysers * dump->nthreads);
  for (int i = 0; i < dump->nthreads; i++) {
    dump->analysers[i] = createOverlapAnalyser(clx, cly);
  }
  return (Dump*) dump;
}
