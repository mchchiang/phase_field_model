// dump.c

#include "phase_field_model.h"
#include "cell.h"
#include "dump.h"

void setDump(Dump* dump, void* derived, char* filename, int printInc,
	     DumpFuncs* funcs) {
  dump->derived = derived;
  dump->filename = filename;
  dump->printInc = printInc;
  dump->funcs = funcs;
}

void deleteDump(Dump* dump) {
  dump->funcs->delete(dump);
  free(dump);
}

void dumpOutput(Dump* dump, PhaseFieldModel* model, int step) {
  if (step % dump->printInc == 0) {
    dump->funcs->output(dump, model, step);
  }
}
