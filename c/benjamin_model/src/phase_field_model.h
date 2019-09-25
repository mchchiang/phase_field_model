// phase_field_model.h

#ifndef PHASE_FIELD_MODEL_H
#define PHASE_FIELD_MODEL_H

#include "cell.h"
#include "dump.h"

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
  double* cellXCM;
  double* cellYCM;
  int* cellXBoundCount;
  int* cellYBoundCount;
  Cell** cells;
  double** totalField;
  Dump** dumps;
  int ndumps;
} PhaseFieldModel;

PhaseFieldModel* createModel(int lx, int ly,int ncells);
void deleteModel(PhaseFieldModel* model);

void initSquareCell(PhaseFieldModel* model, int index,
		    int x, int y, int dx, int dy);

void initCellsFromFile(PhaseFieldModel* model, char* cmFile, char* shapeFile,
		       unsigned long seed);

void run(PhaseFieldModel* model, int nsteps);
void output(PhaseFieldModel* model, int step);

void updateCellField(PhaseFieldModel* model, Cell* cell);
void updateCellVolume(PhaseFieldModel* model, Cell* cell);
void updateCellCM(PhaseFieldModel* model, Cell* cell, int cellIndex);

int iwrap(PhaseFieldModel* model, int i);
int jwrap(PhaseFieldModel* model, int j);
int iup(PhaseFieldModel* model, int i);
int idown(PhaseFieldModel* model, int i);
int jup(PhaseFieldModel* model, int j);
int jdown(PhaseFieldModel* model, int j);
int idiff(PhaseFieldModel* model, int i1, int i2);
int jdiff(PhaseFieldModel* model, int j1, int j2);
int icellwrap(PhaseFieldModel* model, int i);
int jcellwrap(PhaseFieldModel* model, int j);
int icellup(PhaseFieldModel* model, int i);
int icelldown(PhaseFieldModel* model, int i);
int jcellup(PhaseFieldModel* model, int j);
int jcelldown(PhaseFieldModel* model, int j);
double centralDiff(PhaseFieldModel* model, int i, int j, int iu, int id,
		   int ju, int jd, double** field);
double gradient(PhaseFieldModel* model, int i, int j, int u, int d,
		int comp, double** field);
double upwind(PhaseFieldModel* model, int i, int j, int uu, int u,
	      int d, int dd, int comp, double v, double** field);

#endif
