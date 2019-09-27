// dump_energy.c
// Dump the cm of the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "array.h"
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"

typedef struct EnergyDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
  int lx;
  int ly;
  double** phi2Field; // Store the sum of phi^2 field
  double** phi4Field; // Store the sum of phi^4 field
} EnergyDump;

void energyOutput(EnergyDump* dump, PhaseFieldModel* model, int step) {
  FILE* f;
  f = fopen(dump->super.filename, "a");
  
  // Reset fields
  for (int i = 0; i < dump->lx; i++) {
    for (int j = 0; j < dump->ly; j++) {
      dump->phi2Field[i][j] = 0.0;
      dump->phi4Field[i][j] = 0.0;
    }
  }

  // Compute the auxillary field (sum of phi^2)
  Cell* cell;
  double phi, phi2;
  int clx, cly, x, y, cx, cy, get;
  for (int n = 0; n < model->numOfCells; n++) {
    cell = model->cells[n];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    get = cell->getIndex;
    for (int i = 0; i < clx; i++) {
      for (int j = 0; j < cly; j++) {
	x = iwrap(model, cx+i);
	y = jwrap(model, cy+j);
	phi = cell->field[get][i][j];
	phi2 = phi*phi;
	dump->phi2Field[x][y] += phi2;
	dump->phi4Field[x][y] += phi2*phi2;
      }
    }
  }

  // Compute the free energy
  int iu, id, ju, jd;
  double dphi, dphi2, phi2Sum;
  double gradx, grady;
  double cellVolume, dV;
  double energy = 0.0;
  double alphaTerm = 0.0;
  double kappaTerm = 0.0;
  double volumeConst = 0.0;
  double repulsion = 0.0;
  double** cellField;
  for (int n = 0; n < model->numOfCells; n++) {
    cell = model->cells[n];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    get = cell->getIndex;
    cellField = cell->field[get];
    cellVolume = 0.0;
    for (int i = 1; i < clx-1; i++) {
      for (int j = 1; j < cly-1; j++) {
	phi = cellField[i][j];
	phi2 = phi*phi;
	dphi = phi-model->phi0;
	dphi2 = dphi*dphi;
	iu = icellup(model, i);
	id = icelldown(model, i);
	ju = jcellup(model, j);
	jd = jcelldown(model, j);

	// Cahn-Hilliard term
	// Use midpoint method to estimate gradient
	gradx = gradient(model, i, j, iu, id, 0, cellField);
	grady = gradient(model, i, j, ju, jd, 1, cellField);
	alphaTerm += 0.25*model->alpha*phi2*dphi2;
	kappaTerm += 0.5*model->kappa*(gradx*gradx+grady*grady);
	//	cahnHilliard += (0.25*model->alpha*phi2*dphi2 + 
	//		 0.5*model->kappa*(gradx*gradx+grady*grady));
	
	// Calculate cell volume
	cellVolume += phi2;
	
	// Repulsion term
	x = iwrap(model, cx+i);
	y = jwrap(model, cy+j);
	phi2Sum = dump->phi2Field[x][y];
	repulsion += 0.5*model->epsilon*(phi2Sum*phi2Sum-
					 dump->phi4Field[x][y]);
      }
    }
    dV = 1.0-cellVolume/model->piR2phi02;
    volumeConst += model->mu*dV*dV;
  }
  //  energy = cahnHilliard + volumeConst + repulsion;
  energy = alphaTerm + kappaTerm + volumeConst + repulsion;
  fprintf(f, "%d %.5f %.5f %.5f %.5f %.5f\n", step, alphaTerm, kappaTerm,
	  volumeConst, repulsion, energy);
  fclose(f);
}

void deleteEnergyDump(EnergyDump* dump) {
  free(dump->phi2Field);
  free(dump->phi4Field);
  free(dump);
}

DumpFuncs energyDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &energyOutput,
   .destroy = (void (*)(Dump*)) &deleteEnergyDump
  };

Dump* createEnergyDump(char* filename, int lx, int ly, int printInc, 
		       bool overwrite) {
  EnergyDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &energyDumpFuncs);
  dump->overwrite = overwrite;
  dump->lx = lx;
  dump->ly = ly;
  dump->phi2Field = create2DDoubleArray(lx, ly);
  dump->phi4Field = create2DDoubleArray(lx, ly);
  return (Dump*) dump;
}
