// phase_field_model.h

#ifndef PHASE_FIELD_MODEL_H
#define PHASE_FIELD_MODEL_H

#include "cell.h"

typedef struct PhaseFieldModel {
  int lx;
  int ly;
  int numOfCells;
  double dt;
  double piR2phi02;
  double phi0;
  double kappa;
  double alpha;
  double mu;
  double Dr;
  double M;
  double epsilon;
  double motility;
  int cellLx;
  int cellLy;
  Cell** cells;
  double** cellTypeField;
} PhaseFieldModel;

void initModel(PhaseFieldModel* model, int lx, int ly,int ncells);
void deleteModel(PhaseFieldModel* model);

void initSquareCell(PhaseFieldModel* model, int index,
		    int x, int y, int dx, int dy);		    

void run(PhaseFieldModel* model, int nsteps);
void output(PhaseFieldModel* model, int step);

void updateCellField(PhaseFieldModel* model, Cell* cell);
void updateCellVolume(PhaseFieldModel* model);

double singleCellInteractions(PhaseFieldModel* model, Cell* cell,
			      int i, int j);
double cellCellInteractions(PhaseFieldModel* model, Cell* cell, int i, int j);

int iwrap(PhaseFieldModel* model, int i);
int jwrap(PhaseFieldModel* model, int j);
int iup(PhaseFieldModel* model, int i);
int idown(PhaseFieldModel* model, int i);
int jup(PhaseFieldModel* model, int j);
int jdown(PhaseFieldModel* model, int j);
double centralDiff(PhaseFieldModel* model, int i, int j, double** field);
double gradient(PhaseFieldModel* model, int i, int j,
		int comp, double** field);

#endif
