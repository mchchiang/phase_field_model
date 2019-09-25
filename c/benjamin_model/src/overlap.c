// overlap.c

#include <math.h>
#include "array.h"
#include "overlap.h"
#include "cell.h"
#include "phase_field_model.h"

OverlapAnalyser* createOverlapAnalyser(int clx, int cly) {
  OverlapAnalyser* ana = malloc(sizeof *ana);
  ana->cellLx = clx;
  ana->cellLy = cly;
  ana->overlapField = create2DDoubleArray(ana->cellLx, ana->cellLy);
  return ana;
}

void deleteOverlapAnalyser(OverlapAnalyser* ana) {
  free(ana->overlapField);
  free(ana);
}

double computeOverlap(OverlapAnalyser* ana, PhaseFieldModel* model, 
		      int cellIndex) {
  int lx = model->lx;
  int ly = model->ly;
  int clx, cly;
  double olapAvg, product;
  int getm, getn;
  int x0m, y0m, x0n, y0n, dx0mn, dy0mn, xn, yn;
  Cell* cellm;
  Cell* celln;

  // Reset fields to zero
  for (int i = 0; i < ana->cellLx; i++) {
    for (int j = 0; j < ana->cellLy; j++) {
      ana->overlapField[i][j] = 0.0;
    }
  }

  // Compute the sum of all products of two fields: 
  // xi_i(x) = sum_{j}phi_i(x)phi_j(x)
  olapAvg = 0.0;
  cellm = model->cells[cellIndex];
  getm = cellm->getIndex;
  clx = cellm->lx;
  cly = cellm->ly;
  x0m = iwrap(model, cellm->x);
  y0m = jwrap(model, cellm->y);
  for (int n = 0; n < model->numOfCells; n++) {
    celln = model->cells[n];
    getn = celln->getIndex;
    if (cellIndex == n) continue;
    x0n = iwrap(model, celln->x);
    y0n = jwrap(model, celln->y);
    dx0mn = idiff(model, x0m, x0n);
    dy0mn = jdiff(model, y0m, y0n);
    if (abs(dx0mn) > lx || abs(dy0mn) > ly) continue;
    for (int i = 0; i < clx; i++) {
      xn = dx0mn+i;
      if (xn < 0 || xn >= clx) continue;
      for (int j = 0; j < cly; j++) {
	yn = dy0mn+j;
	if (yn < 0 || yn >= cly) continue;
	product = cellm->field[getm][i][j]*celln->field[getn][xn][yn];
	ana->overlapField[i][j] += product;
	olapAvg += product;
      }
    }
  }
  olapAvg /= (lx*ly);
  return olapAvg;
}
