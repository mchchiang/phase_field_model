// test_shape.c

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "array.h"
#include "shape.h"

int main(int argc, char* argv[]) {
  
  int n = atoi(argv[1]);
  double** field = create2DDoubleArray(n, n);
  ShapeAnalyser* ana = createShapeAnalyser(4, n, n, 25, 4.0, 60, 101, 1.0);
  
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
  getShapeInfo(ana, field, &perimeter, &area, &shapeIndex);
  printf("perimeter = %.5f area = %.5f ratio = %.5f\n",
	 perimeter, area, shapeIndex);
  
  free(field);
  deleteShapeAnalyser(ana);
}
