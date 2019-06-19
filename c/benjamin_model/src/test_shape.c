// test_shape.c

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "array.h"
#include "shape.h"

int main(int argc, char* argv[]) {
  
  int n = atoi(argv[1]);
  double** field = create2DDoubleArray(n, n);
  int nthreads = omp_get_max_threads();
  ShapeAnalyser** ana = malloc(sizeof *ana * nthreads);
  for (int i = 0; i < nthreads; i++) {
    ana[i] = createShapeAnalyser(4, n, n, 25, 4.0, 60, 101, 1.0);
  }
  
  // Read the shape file
  char line [80];
  FILE* f = fopen(argv[2], "r");
  int x, y;
  double value;
  while (fgets(line, sizeof(line), f) != NULL) {
    if (sscanf(line, "%d %d %lf", &x, &y, &value)) {
      field[x][y] = value;
    }
  }
  fclose(f);
  
  double perimeter, area, shapeIndex;
  int i;
  int id = 0;
#pragma omp parallel default(none) shared(ana, field) \
  private (i, id, perimeter, area, shapeIndex)
{
  id = omp_get_thread_num();
#pragma omp for schedule(static)
  for (i = 0; i < 100; i++) {
    getShapeInfo(ana[id], field, &perimeter, &area, &shapeIndex);
    printf("Done %d: perimeter = %.5f area = %.5f ratio = %.5f\n",
	   i, perimeter, area, shapeIndex);
  }
}
  free(field);
  for (i = 0; i < nthreads; i++) {
    deleteShapeAnalyser(ana[i]);
  }
  free(ana);
}
