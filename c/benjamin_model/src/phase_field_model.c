// phase_field_model.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "array.h"
#include "cell.h"
#include "phase_field_model.h"
#include "dump.h"

PhaseFieldModel* createModel(int _lx, int _ly, int ncells) {
  PhaseFieldModel* model = (PhaseFieldModel*) malloc(sizeof(PhaseFieldModel));
  model->lx = _lx;
  model->ly = _ly;
  model->numOfCells = ncells;
  model->dt = 0.01;
  model->phi0 = 1.0;
  model->piR2phi02 = M_PI;
  model->kappa = 0.0;
  model->alpha = 0.0;
  model->mu = 0.0;
  model->Dr = 1.0;
  model->M = 1.0;
  model->epsilon = 0.0;
  model->motility = 0.0;
  model->cellLx = 1;
  model->cellLy = 1;
  model->cells = (Cell**) malloc(sizeof(Cell*)*ncells);
  model->cellTypeField = create2DDoubleArray(model->lx, model->ly);
  model->dump[0] = createFieldDump("field.dat", 1000);
  model->dump[1] = createCMDump("cm.dat", 1000);
  return model; 
}

void deleteModel(PhaseFieldModel* model) {
  for (int i = 0; i < model->numOfCells; i++) {
    if (model->cells[i] != NULL) {
      deleteCell(model->cells[i]);
    }
  }
  free(model->cells);
  free(model->cellTypeField);
  for (int i = 0; i < 2; i++) {
    deleteDump(model->dump[i]);
  }
  free(model);
}

void initSquareCell(PhaseFieldModel* model, int index, int x, int y,
		    int dx, int dy) {
  int clx = model->cellLx;
  int cly = model->cellLy;
  Cell* cell = initCell(x, y, clx, cly, model->phi0/2.0, index);
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
      cell = initCell(x, y, clx, cly, model->phi0/2.0, index+seed);
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
  printf("Running model ...\n");
  for (int n = 1; n <= nsteps; n++) {
    // Update each cell type field
    // Reset fields to zero
    for (int i = 0; i < model->lx; i++) {
      for (int j = 0; j < model->ly; j++) {
	model->cellTypeField[i][j] = 0.0;
      }
    }
    
    Cell* cell;
    double phi; 
    int clx, cly, x, y, cx, cy, get;
    for (int i = 0; i < model->numOfCells; i++) {
      cell = model->cells[i];
      clx = cell->lx;
      cly = cell->ly;
      cx = cell->x;
      cy = cell->y;
      get = cell->getIndex;
      for (int j = 0; j < clx; j++) {
	for (int k = 0; k < cly; k++) {
	  x = iwrap(model, cx+j);
	  y = jwrap(model, cy+k);
	  phi = cell->field[get][j][k];
	  model->cellTypeField[x][y] += phi*phi;
	}
      }
    }
    
    // Update each cell field
    int i;
#pragma omp parallel default(none) shared(model) private(i, cell) 
    {
#pragma omp for schedule(dynamic, 4)
    for (i = 0; i < model->numOfCells; i++) {
      cell = model->cells[i];
      updateVolume(cell);
      updateCellField(model, cell);
      updateCM(cell);
      updateVelocity(cell);
    }
    }
    output(model, n);
  }
}

void output(PhaseFieldModel* model, int step) {
  if (step % 1000 == 0) {
    printf("Step %d\n", step);
  }
  for (int i = 0; i < 2; i++) {
    dumpOutput(model->dump[i], model, step);
  }
}

void updateCellField(PhaseFieldModel* model, Cell* cell) {
  startUpdateCellField(cell);
  int cx = cell->lx;
  int cy = cell->ly;
  // Apply fixed (Dirichlet) boundary conditions (u = 0 at boundaries)
  // i and j are coordinates in the cell's own reference frame
  int set = cell->setIndex;
  int get = cell->getIndex;
  for (int i = 1; i < cx-1; i++) {
    for (int j = 1; j < cy-1; j++) {
      cell->field[set][i][j] = cell->field[get][i][j] + model->dt * model->M *
	(singleCellInteractions(model, cell, i, j) +
	 cellCellInteractions(model, cell, i, j));
    }
  }
  endUpdateCellField(cell);
}

double singleCellInteractions(PhaseFieldModel* model, Cell* cell,
			      int i, int j) {
  double** cellField = cell->field[cell->getIndex];
  double phi = cellField[i][j];
  double cahnHilliard = model->kappa * centralDiff(model, i, j, cellField) +
    model->alpha * phi * (model->phi0 - phi) * (phi - 0.5 * model->phi0);
  double advection =  -model->motility / model->M * 
    (cell->vx * gradient(model, i, j, 0, cellField) -
     cell->vy * gradient(model, i, j, 1, cellField));
  double fixVolume = 4.0 * model->mu * phi / (model->piR2phi02) *
    (1.0 - 1.0/(model->piR2phi02) * cell->volume);
  return cahnHilliard + fixVolume + advection;
}

double cellCellInteractions(PhaseFieldModel* model, Cell* cell, int i, int j) {
  double phi = cell->field[cell->getIndex][i][j];
  int cx = cell->x;
  int cy = cell->y;
  int x = iwrap(model, cx+i);
  int y = jwrap(model, cy+j);

  // Compute exclusion effect
  return 2.0 * model->epsilon * phi * (phi*phi - model->cellTypeField[x][y]);
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

int iup(PhaseFieldModel* model, int i) {
  return (i+1 >= model->lx) ? 0 : i+1;
}

int idown(PhaseFieldModel* model, int i) {
  return (i-1 < 0) ? model->lx-1 : i-1;
}

int jup(PhaseFieldModel* model, int j) {
  return (j+1 >= model->ly) ? 0 : j+1;
}

int jdown(PhaseFieldModel* model, int j) {
  return (j-1 < 0) ? model->ly-1 : j-1;
}

double centralDiff(PhaseFieldModel* model, int i, int j, double** field) {
  return field[iup(model,i)][j] + field[idown(model,i)][j] +
    field[i][jup(model,j)] + field[i][jdown(model,j)] - 4.0 * field[i][j];
}

double gradient(PhaseFieldModel* model, int i, int j,
		int comp, double** field) {
  if (comp == 0) {
    return 0.5*(field[iup(model,i)][j] - field[idown(model,i)][j]);
  } else {
    return 0.5*(field[i][jup(model,j)] - field[i][jdown(model,j)]);
  }
}
