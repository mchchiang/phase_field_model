// phase_field_model.c

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include "array.h"
#include "cell.h"
#include "phase_field_model.h"
#include "dump.h"

PhaseFieldModel* createModel(int lx, int ly, int ncells) {
  PhaseFieldModel* model =  malloc(sizeof *model);
  model->lx = lx;
  model->ly = ly;
  model->numOfCells = ncells;
  model->dt = 0.01;
  model->phi0 = 1.0;
  model->piR2phi02 = M_PI;
  model->kappa = 0.0;
  model->alpha = 0.0;
  model->mu = 0.0;
  model->Dr = 0.0;
  model->M = 1.0;
  model->epsilon = 0.0;
  model->motility = 0.0;
  model->cellLx = 1;
  model->cellLy = 1;
  model->cells = malloc(sizeof *model->cells * ncells);
  model->totalField = create2DDoubleArray(model->lx, model->ly);
  model->dumps = NULL;
  model->ndumps = 0;
  return model; 
}

void deleteModel(PhaseFieldModel* model) {
  for (int i = 0; i < model->numOfCells; i++) {
    if (model->cells[i] != NULL) {
      deleteCell(model->cells[i]);
    }
  }
  free(model->cells);
  free(model->totalField);
  free(model);
}

void initSquareCell(PhaseFieldModel* model, int index, int x, int y,
		    int dx, int dy) {
  int clx = model->cellLx;
  int cly = model->cellLy;
  Cell* cell = createCell(x, y, clx, cly, model->Dr, model->phi0/2.0, index);
  model->cells[index] = cell;
  int x0 = (clx-dx)/2;
  int y0 = (cly-dy)/2;
  initFieldSquare(cell, x0, y0, dx, dy, model->phi0);
}

void initCellsFromFile(PhaseFieldModel* model, char* cmFile,
		       char* shapeFile, int seed) {
  FILE* fcm = fopen(cmFile, "r");
  if (fcm == NULL) {
    printf("Problem with opening the centre of mass file!\n");
    return;
  }

  FILE* fshape = fopen(shapeFile, "r");
  if (fshape == NULL) {
    printf("Problem with opening the shape file!\n");
    return;
  }

  char line [80];
  int nvar, x, y;
  double val;
  int clx = model->cellLx;
  int cly = model->cellLy;

  // Allocate memory for the template field
  double** field = create2DDoubleArray(clx, cly);
  // Read the template field from the shape file
  while (fgets(line, sizeof(line), fshape) != NULL) {
    nvar = sscanf(line, "%d %d %lf", &x, &y, &val);
    if (nvar == 3 && x >= 0 && x < clx && y >= 0 && y < cly) {
      field[x][y] = val;
    }
  }
  
  int index;
  double xcm, ycm;
  int count = 0;
  Cell* cell;
  while (fgets(line, sizeof(line), fcm) != NULL) {
    nvar = sscanf(line, "%d %lf %lf", &index, &xcm, &ycm);
    if (nvar == 3) {
      x = (int)(xcm-clx/2);
      y = (int)(ycm-cly/2);
      cell = createCell(x, y, clx, cly, model->Dr, 
			model->phi0/2.0, index+seed);
      model->cells[index] = cell;
      initField(cell, field);
      count++;
    }
  }

  if (count != model->numOfCells) {
    printf("Not all cells initialised!\n");
  }
  free(field);
  fclose(fcm);
  fclose(fshape);
}

void run(PhaseFieldModel* model, int nsteps) {
  output(model, 0); // Output initial state of the model
  for (int n = 1; n <= nsteps; n++) {
    // Update total field
    for (int i = 0; i < model->lx; i++) {
      for (int j = 0; j < model->ly; j++) {
	model->totalField[i][j] = 0.0;
      }
    }
    
    Cell* cell;
    double phi; 
    int clx, cly, x, y, cx, cy, get;
    int i, j, k;
#pragma omp parallel default(none) shared(model) \
private(i, j, k, clx, cly, x, y, cx, cy, get, phi, cell) 
{
#pragma omp for schedule(static)
    for (i = 0; i < model->numOfCells; i++) {
      cell = model->cells[i];
      clx = cell->lx;
      cly = cell->ly;
      cx = cell->x;
      cy = cell->y;
      get = cell->getIndex;
      for (j = 0; j < clx; j++) {
	for (k = 0; k < cly; k++) {
	  x = iwrap(model, cx+j);
	  y = jwrap(model, cy+k);
	  phi = cell->field[get][j][k];
#pragma omp atomic
	  model->totalField[x][y] += phi*phi;
	}
      }
    }
    // Implicit sync barrier after for loop

    // Update each cell field    
#pragma omp for schedule(static)
    for (i = 0; i < model->numOfCells; i++) {
      cell = model->cells[i];
      updateVolume(cell);
      updateCellField(model, cell);
      updateCM(cell);
      updateGyration(cell);
      updateVelocity(cell, model->dt);
    }
    }
    output(model, n);
  }
}

void output(PhaseFieldModel* model, int step) {
  if (step % 1000 == 0) {
    printf("Step %d\n", step);
  }
  for (int i = 0; i < model->ndumps; i++) {
    dumpOutput(model->dumps[i], model, step);
  }
}

void updateCellField(PhaseFieldModel* model, Cell* cell) {
  startUpdateCellField(cell);
  int clx = cell->lx;
  int cly = cell->ly;
  int cx = cell->x;
  int cy = cell->y;
  int set = cell->setIndex;
  int get = cell->getIndex;

  double** cellField = cell->field[get];
  double cahnHilliard, volumeConst, advection, repulsion, phi;
  int x, y; // Lab frame coordinates of a lattice element
  int iuu, iu, id, idd, juu, ju, jd, jdd;  // Nearest neighbours
  double vol = cell->volume / model->piR2phi02;
  double volprefactor = 4.0 * model->mu / model->piR2phi02;

  // Apply fixed (Dirichlet) boundary conditions (u = 0 at boundaries)
  // i and j are coordinates in the cell's own reference frame
  for (int i = 2; i < clx-2; i++) {
    for (int j = 2; j < cly-2; j++) {
      phi = cellField[i][j];
      iu = icellup(model, i);
      iuu = icellup(model, iu);
      id = icelldown(model, i);
      idd = icelldown(model, id);
      ju = jcellup(model, j);
      juu = jcellup(model, ju);
      jd = jcelldown(model, j);
      jdd = jcelldown(model, jd);
      
      // Cahn-Hilliard term
      cahnHilliard = model->kappa *
	centralDiff(model, i, j, iu, id, ju, jd, cellField) +
	model->alpha * phi * (model->phi0 - phi) * (phi - 0.5 * model->phi0);

      // Volume term
      volumeConst = volprefactor * phi * (1.0 - vol);

      // Advection term (use the 3rd order upwind scheme)
      advection = model->motility *
	(upwind(model, i, j, iuu, iu, id, idd, 0, cell->vx, cellField) +
	 upwind(model, i, j, juu, ju, jd, jdd, 1, cell->vy, cellField));

      // Repulsion term
      x = iwrap(model, cx+i);
      y = jwrap(model, cy+j);
      repulsion = 2.0 * model->epsilon * phi *
	(phi*phi - model->totalField[x][y]);
      
      cell->field[set][i][j] = cell->field[get][i][j] + model->dt *
	(model->M * (cahnHilliard + volumeConst + repulsion) - advection);
    }
  }
  endUpdateCellField(cell);
}

int iwrap(PhaseFieldModel* model, int i) {
  int remainder = i % model->lx;
  if (remainder >= 0) {
    return remainder;
  }
  return model->lx + remainder;
}

int jwrap(PhaseFieldModel* model, int j) {
  int remainder = j % model->ly;
  if (remainder >= 0) {
    return remainder;
  }
  return model->ly + remainder;
}

inline int iup(PhaseFieldModel* model, int i) {
  return (i+1 >= model->lx) ? 0 : i+1;
}

inline int idown(PhaseFieldModel* model, int i) {
  return (i-1 < 0) ? model->lx-1 : i-1;
}

inline int jup(PhaseFieldModel* model, int j) {
  return (j+1 >= model->ly) ? 0 : j+1;
}

inline int jdown(PhaseFieldModel* model, int j) {
  return (j-1 < 0) ? model->ly-1 : j-1;
}

int icellwrap(PhaseFieldModel* model, int i) {
  int remainder = i % model->cellLx;
  if (remainder >= 0) {
    return remainder;
  }
  return model->cellLx + remainder;
}

int jcellwrap(PhaseFieldModel* model, int j) {
  int remainder = j % model->cellLy;
  if (remainder >= 0) {
    return remainder;
  }
  return model->cellLy + remainder;
}

inline int icellup(PhaseFieldModel* model, int i) {
  return (i+1 >= model->cellLx) ? 0 : i+1;
}

inline int icelldown(PhaseFieldModel* model, int i) {
  return (i-1 < 0) ? model->cellLx-1 : i-1;
}

inline int jcellup(PhaseFieldModel* model, int j) {
  return (j+1 >= model->cellLy) ? 0 : j+1;
}

inline int jcelldown(PhaseFieldModel* model, int j) {
  return (j-1 < 0) ? model->cellLy-1 : j-1;
}

inline double centralDiff(PhaseFieldModel* model, int i, int j, int iu, int id,
		   int ju, int jd, double** field) {
  return field[iu][j] + field[id][j] + field[i][ju] + field[i][jd] -
    4.0 * field[i][j];
}

inline double gradient(PhaseFieldModel* model, int i, int j, int u, int d,
		       int comp, double** field) {
  switch (comp) {
  case 0: return 0.5*(field[u][j] - field[d][j]);
  case 1: return 0.5*(field[i][u] - field[i][d]);
  default: return 0.0;
  }
}

inline double upwind(PhaseFieldModel* model, int i, int j, int uu, int u,
		     int d, int dd, int comp, double v, double** field) {
  switch (comp) {
  case 0: return v > 0.0 ?
      v * (2.0 * field[u][j] + 3.0 * field[i][j] -
	   6.0 * field[d][j] + field[dd][j]) / 6.0 :
    v * (-field[uu][j] + 6.0 * field[u][j] -
	 3.0 * field[i][j] - 2.0 * field[d][j]) / 6.0;
  case 1: return v > 0.0 ?
      v * (2.0 * field[i][u] + 3.0 * field[i][j] -
	   6.0 * field[i][d] + field[i][dd]) / 6.0 :
    v * (-field[i][uu] + 6.0 * field[i][u] -
	 3.0 * field[i][j] - 2.0 * field[i][d]) / 6.0;
  default: return 0.0;
  }
}
