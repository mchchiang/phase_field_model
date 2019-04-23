// dump.c

#include <string.h>
#include "phase_field_model.h"
#include "cell.h"
#include "dump.h"

void setDump(Dump* dump, void* derived, char* filename, int printInc,
	     DumpFuncs* funcs) {
  dump->derived = derived;
  char* filenamecpy = malloc(sizeof *filenamecpy * DIR_SIZE);
  strcpy(filenamecpy, filename);
  dump->filename = filenamecpy;
  dump->printInc = printInc;
  dump->funcs = funcs;
}

void deleteDump(Dump* dump) {
  dump->funcs->delete(dump);
  free(dump->filename);
  free(dump);
}

void dumpOutput(Dump* dump, PhaseFieldModel* model, int step) {
  if (step % dump->printInc == 0) {
    dump->funcs->output(dump, model, step);
  }
}
