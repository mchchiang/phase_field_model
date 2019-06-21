// dump.h

#ifndef DUMP_H
#define DUMP_H

#define DIR_SIZE 500

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

struct PhaseFieldModel;
struct Dump;

typedef struct DumpFuncs {
  void (*output)(struct Dump* dump, struct PhaseFieldModel* model, int step);
  void (*destroy)(struct Dump* dump);
} DumpFuncs;

typedef struct Dump {
  DumpFuncs* funcs;
  void* derived;
  char* filename;
  int printInc;
} Dump;

void dumpOutput(Dump* dump, struct PhaseFieldModel* model, int step);
void deleteDump(Dump* dump);
void setDump(Dump* dump, void* derived, char* filename, int printInc,
	     DumpFuncs* funcs);

Dump* createCellFieldDump(char* filename, int cellIndex,
			  int printInc, bool overwrite);
Dump* createFieldDump(char* filename, int printInc, bool overwrite);
Dump* createCMDump(char* filename, int printInc, bool overwrite);
Dump* createBulkCMDump(char* filenmae, int printInc);
Dump* createGyrationDump(char* filename, int printInc, bool overwrite);
Dump* createNeighbourDump(char* filename, int lx, int ly, int printInc,
			  bool overwrite);
Dump* createShapeDump(char* filename, int scale, int lx, int ly,
		      int kernelLength, double sigma, int sgolayDegree,
		      int sgolayLength, double threshold, int printInc,
		      bool overwrite);

#endif
