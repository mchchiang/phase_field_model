// phase_field_model.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "array.h"
#include "cell.h"
#include "phase_field_model.h"

void initModel(PhaseFieldModel* model, int _lx, int _ly,
	       int _ntypes, int _ncells) {
  model->lx = _lx;
  model->ly = _ly;
  int ntypes = _ntypes;
  int ncells = _ncells;
  model->dt = 0.01;
  model->numOfCellTypes = ntypes;
  model->numOfCells = ncells;
  model->cells = (Cell**) malloc(sizeof(Cell*)*ncells);
  for (int i = 0; i < ncells; i++) {
    model->cells[i] = (Cell*) malloc(sizeof(Cell));
  }
  model->cellTypeField = (double***) malloc(sizeof(double**)*ntypes);
  for (int i = 0; i < ntypes; i++) {
    model->cellTypeField[i] = create2DDoubleArray(model->lx, model->ly);
  }
  model->idealCellVolume = create1DDoubleArray(ntypes);
  model->diffusion = create1DDoubleArray(ntypes);
  int ntypesSq = ntypes*ntypes;
  model->alpha = create1DDoubleArray(ntypes);
  model->beta = create1DDoubleArray(ntypesSq);
  model->eta = create1DDoubleArray(ntypesSq);
  model->gamma = create1DDoubleArray(ntypes);
  model->motility = create1DDoubleArray(ntypes);
  model->cellLx = create1DIntArray(ntypes);
  model->cellLy = create1DIntArray(ntypes);
}

void deleteModel(PhaseFieldModel* model) {
  for (int i = 0; i < model->numOfCells; i++) {
    if (model->cells[i] != NULL) {
      deleteCell(model->cells[i]);
    }
  }
  free(model->cells);
  for (int i = 0; i < model->numOfCellTypes; i++) {
    free(model->cellTypeField[i]);
  }
  free(model->cellTypeField);
  free(model->idealCellVolume);
  free(model->alpha);
  free(model->beta);
  free(model->eta);
  free(model->gamma);
  free(model->cellLx);
  free(model->cellLy);
  free(model->motility);
}

void initSquareCell(PhaseFieldModel* model, int index, int x, int y,
		    int dx, int dy, int type) {
  Cell* cell = model->cells[index];
  int clx = model->cellLx[type];
  int cly = model->cellLy[type];
  initCell(cell, x, y, clx, cly, type, index);
  int x0 = (clx-dx)/2;
  int y0 = (cly-dy)/2;
  initFieldSquareOnes(cell, x0, y0, dx, dy);
}

void run(PhaseFieldModel* model, int nsteps) {
  printf("Running model ...\n");
  for (int n = 0; n < nsteps; n++) {
    // Update each cell type field
    for (int i = 0; i < model->numOfCellTypes; i++) {
      // Reset fields to zero
      for (int j = 0; j < model->lx; j++) {
	for (int k = 0; k < model->ly; k++) {
	  model->cellTypeField[i][j][k] = 0.0;
	}
      }
    }
    
    Cell* cell;
    int t, clx, cly, x, y, cx, cy;
    for (int i = 0; i < model->numOfCells; i++) {
      cell = model->cells[i];
      t = cell->type;
      clx = cell->lx;
      cly = cell->ly;
      cx = cell->x;
      cy = cell->y;
      for (int j = 0; j < clx; j++) {
	for (int k = 0; k < cly; k++) {
	  x = iwrap(model, cx+j);
	  y = jwrap(model, cy+k);
	  model->cellTypeField[t][x][y] += cell->volumeField[j][k];
	}
      }
    }
    
    // Update each cell field
    //    int i, j, k, set, get;
    // double u, double totalVolume;
    int i;
#pragma omp parallel default(none) shared(model) private(i, cell) 
    {
#pragma omp for schedule(dynamic, 1)
    for (i = 0; i < model->numOfCells; i++) {
      cell = model->cells[i];
      updateVolume(cell);
      updateCellField(model, cell);
      updateCM(cell);
      updateVelocity(cell);
    }
    }
    output(model, n, 0);
  }
}

void output(PhaseFieldModel* model, int step, int type) {
  if (step % 1000 == 0) {
    printf("Step %d\n", step);
    FILE* f = fopen("output.dat.tmp", "w");
    for (int i = 0; i < model->lx; i++) {
      for (int j = 0; j < model->ly; j++) {
	fprintf(f, "%.5f ", model->cellTypeField[type][i][j]);
      }
      fprintf(f, "\n");
    }
    fclose(f);
    //    char buf[80];
    //sprintf(buf, "output.dat.%d", step);
    rename("output.dat.tmp", "output.dat");
    
    FILE* fcm = fopen("cm.dat.tmp", "w");
    Cell* cell;
    double x, y, cx, cy;
    for (int i = 0; i < model->numOfCells; i++) {
      cell = model->cells[i];
      cx = cell->x;
      cy = cell->y;
      x = iwrap(model, cx+cell->xcm);
      y = jwrap(model, cy+cell->ycm);
      fprintf(fcm, "%d %.5f %.5f\n", i, x, y);
    }
    fclose(fcm);
    rename("cm.dat.tmp", "cm.dat");
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
      cell->field[set][i][j] = cell->field[get][i][j] + model->dt *
	(singleCellInteractions(model, cell, i, j) +
	 cellCellInteractions(model, cell, i, j) +
	 cellSubstrateInteractions(model, cell, i, j));
    }
  }
  endUpdateCellField(cell);
}

double singleCellInteractions(PhaseFieldModel* model, Cell* cell,
			      int i, int j) {
  double** cellField = cell->field[cell->getIndex];
  double u = cellField[i][j];
  int t = cell->type;
  return model->diffusion[t] * centralDiff(model, i, j, cellField) +
    -model->motility[t] * (cell->vx * gradient(model, i, j, 0, cellField) -
			   cell->vy * gradient(model, i, j, 1, cellField)) +
    u * (1 - u) * (u - 0.5 + model->alpha[t] *
		   (model->idealCellVolume[t] - cell->volume));
}

double cellCellInteractions(PhaseFieldModel* model, Cell* cell, int i, int j) {
  double u = cell->field[cell->getIndex][i][j];
  int t = cell->type;
  int cx = cell->x;
  int cy = cell->y;
  int x = iwrap(model, cx+i);
  int y = jwrap(model, cy+j);

  // Compute excluded volume effect
  double exclusion = model->beta[getTypeIndex(model, t, t)] *
    cell->volumeField[i][j];
  for (int l = 0; l < model->numOfCellTypes; l++) {
    exclusion -= model->beta[getTypeIndex(model, l, t)] *
      model->cellTypeField[l][x][y];
  }
  
  // Compute adhesion effect
  double volumeFieldCentralDiff = centralDiff(model, i, j, cell->volumeField);
  double adhesion = -model->eta[getTypeIndex(model, t, t)] *
    volumeFieldCentralDiff;
  for (int l = 0; l < model->numOfCellTypes; l++) {
    adhesion += model->eta[getTypeIndex(model, l, t)] *
      centralDiff(model, x, y, model->cellTypeField[l]);
  }

  // Compute regularisation effect
  double regularisation = model->gamma[t] * volumeFieldCentralDiff;
  return u * (1 - u) * (exclusion + adhesion + regularisation);
}

double cellSubstrateInteractions(PhaseFieldModel* model, Cell* cell,
				 int i, int j) {
  return 0.0; // Not implemented at the moment
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

int getTypeIndex(PhaseFieldModel* model, int type1, int type2) {
  return type1*model->numOfCellTypes+type2;
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

double inCellFunc(int i, int j, double** field) {
  double u = field[i][j];
  return u*u*(3-2*u);
}

