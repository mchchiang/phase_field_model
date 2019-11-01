// dump_shape.c
// Dump the shape info (area, perimeter, and asphericity) of the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"
#include "shape.h"
#include "array.h"
#include "constant.h"

typedef struct ShapeDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
  int nthreads; // Use multiple threads to process the images
  ShapeAnalyser** analysers; // One analyser per thread
} ShapeDump;

void shapeOutput(ShapeDump* dump, PhaseFieldModel* model, int step) {
  char tmpfile [PF_DIR_SIZE];
  FILE* f;
  if (dump->overwrite) {
    strcpy(tmpfile, dump->super.filename);
    strcat(tmpfile, ".tmp");
    f = fopen(tmpfile, "w");
  } else {
    f = fopen(dump->super.filename, "a");
  }
  // Compute the area, perimeter and asphericity of each cell
  int j;
  int id = 0;
  int ncells = model->numOfCells;
  ShapeInfo** info = malloc(sizeof *info * ncells);
  Cell* cell;
  ShapeAnalyser* analyser;
#pragma omp parallel default(none) \
  shared(ncells, dump, model, info)	\
  private (j, id, cell, analyser)
  {
#ifdef _OPENMP
    id = omp_get_thread_num();
#endif
    analyser = dump->analysers[id];
#pragma omp for schedule(static)
    for (j = 0; j < ncells; j++) {
      cell = model->cells[j];
      info[j] = getShapeInfo(analyser, cell->field[cell->getIndex]);
    }
  }
  // Output the shape data to file
  fprintf(f, "Cells: %d\n", model->numOfCells);
  fprintf(f, "Timestep: %d\n", step);
  int pixels;
  double perimeter, area, pixelArea, chainPerimeter;
  for (int i = 0; i < ncells; i++) {
    pixels = info[i]->pixels;
    area = info[i]->area;
    perimeter = info[i]->perimeter;
    pixelArea = info[i]->pixelArea;
    chainPerimeter = info[i]->chainPerimeter;
    fprintf(f, "%.5f %.5f %.5f %.5f %d\n", perimeter, area,
	    chainPerimeter, pixelArea, pixels);
  }
  fclose(f);
  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }

  // Clean up resources
  for (int i = 0; i < ncells; i++) {
    deleteShapeInfo(info[i]);
  }
  free(info);
}

void deleteShapeDump(ShapeDump* dump) {
  for (int i = 0; i < dump->nthreads; i++) {
    deleteShapeAnalyser(dump->analysers[i]);
  }
  free(dump->analysers);
  free(dump);
}

DumpFuncs shapeDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &shapeOutput,
   .destroy = (void (*)(Dump*)) &deleteShapeDump
  };

Dump* createShapeDump(char* filename, int scale, int lx, int ly,
		      int kernelLength, double sigma, int sgolayDegree,
		      int sgolayLength, double threshold,
		      int printInc, bool overwrite) {
  ShapeDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &shapeDumpFuncs);
  dump->overwrite = overwrite;
  // Initialise the shape analysers
#ifdef _OPENMP
  dump->nthreads = omp_get_max_threads();
#else
  dump->nthreads = 1;
#endif
  dump->analysers = malloc(sizeof *dump->analysers * dump->nthreads);
  for (int i = 0; i < dump->nthreads; i++) {
    dump->analysers[i] = createShapeAnalyser(scale, lx, ly, kernelLength,
					     sigma, sgolayDegree, sgolayLength,
					     threshold);
  }
  return (Dump*) dump;
}
