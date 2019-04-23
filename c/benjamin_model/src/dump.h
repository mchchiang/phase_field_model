// dump.h

#ifndef DUMP_H
#define DUMP_H

#define DUMP_DIR_SIZE 500
#define DUMP_FILE_SIZE 80

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

struct PhaseFieldModel;
struct Dump;

typedef struct DumpFuncs {
  void (*output)(struct Dump* dump, struct PhaseFieldModel* model, int step);
  void (*delete)(struct Dump* dump);
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
			  int printInc, bool override);
Dump* createFieldDump(char* filename, int printInc, bool override);
Dump* createCMDump(char* filename, int printInc, bool override);
Dump* createGyrationDump(char* filename, int printInc, bool override);

#endif
