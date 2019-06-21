// dump.c

#include <string.h>
#include "phase_field_model.h"
#include "cell.h"
#include "dump.h"

void setDump(Dump* dump, void* derived, char* filename, int printInc,
	     DumpFuncs* funcs) {
  dump->derived = derived;
  dump->filename = malloc(sizeof *dump->filename * PF_DIR_SIZE);
  strcpy(dump->filename, filename);
  dump->printInc = printInc;
  dump->funcs = funcs;
}

void deleteDump(Dump* dump) {
  free(dump->filename);
  dump->funcs->destroy(dump);
}

void dumpOutput(Dump* dump, PhaseFieldModel* model, int step) {
  if (step % dump->printInc == 0) {
    dump->funcs->output(dump, model, step);
  }
}
