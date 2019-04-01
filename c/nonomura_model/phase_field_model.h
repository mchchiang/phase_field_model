// phase_field_model.h

#ifndef PHASE_FIELD_MODEL_H
#define PHASE_FIELD_MODEL_H

#include "cell.h"

typedef struct PhaseFieldModel {
  int lx;
  int ly;
  int numOfCellTypes;
  int numOfCells;
  double dt;
  Cell** cells;
  double*** cellTypeField;
  double* idealCellVolume;
  double* diffusion; 
  double* alpha;
  double* beta;
  double* eta;
  double* gamma;
  double* motility;
  int* cellLx;
  int* cellLy;
} PhaseFieldModel;

void initModel(PhaseFieldModel* model, int lx, int ly, int ntypes, int ncells);
void deleteModel(PhaseFieldModel* model);

void initSquareCell(PhaseFieldModel* model, int index,
		    int x, int y, int dx, int dy, int type);
		    

void run(PhaseFieldModel* model, int nsteps);
void output(PhaseFieldModel* model, int step, int type);

void updateCellField(PhaseFieldModel* model, Cell* cell);
void updateCellVolume(PhaseFieldModel* model);

double singleCellInteractions(PhaseFieldModel* model, Cell* cell,
			      int i, int j);
double cellCellInteractions(PhaseFieldModel* model, Cell* cell, int i, int j);
double cellSubstrateInteractions(PhaseFieldModel* model, Cell* cell,
				 int i, int j);

int iwrap(PhaseFieldModel* model, int i);
int jwrap(PhaseFieldModel* model, int j);
int iup(PhaseFieldModel* model, int i);
int idown(PhaseFieldModel* model, int i);
int jup(PhaseFieldModel* model, int j);
int jdown(PhaseFieldModel* model, int j);
int getTypeIndex(PhaseFieldModel* model, int type1, int type2);
double centralDiff(PhaseFieldModel* model, int i, int j, double** field);
double gradient(PhaseFieldModel* model, int i, int j,
		int comp, double** field);

#endif
