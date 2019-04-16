// dump.h

#ifndef DUMP_H
#define DUMP_H

#include <stdio.h>
#include <stdlib.h>

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

Dump* createCellFieldDump(char* filename, int cellIndex, int printInc);
Dump* createFieldDump(char* filename, int printInc);
Dump* createCMDump(char* filename, int printInc);

#endif
