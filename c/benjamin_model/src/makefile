
CC := gcc-8
# CC := gcc
CFLAGS := -Wall -Wshadow -Winline -O3 -std=c11 -fopenmp
INC := -I../../mtwister/ -I../../arralloc/
LDFLAGS := -lm -lgomp

PRG_SUFFIX_FLAG := 0

PRG_SUFFIX =

ODIR := ../bin/obj/
EDIR := ../bin/exe/

SRC := cell.c phase_field_model.c array.c dump.c dump_cm.c dump_field.c dump_cell_field.c dump_gyration.c random_mtwister.c
OBJ := $(patsubst %.c, $(ODIR)/%.o, $(SRC))
PRG_SRC := run_phase_field_model.c
PRG_OBJ := $(patsubst %.c, $(ODIR)/%.o, $(PRG_SRC))
LSRC := ../../mtwister/mtwister.c ../../arralloc/arralloc.c
LOBJ := $(patsubst %.c, %.o, $(LSRC)) #../mtwister/mtwister.o ../arralloc/arralloc.o
PRG := $(patsubst %.c, $(EDIR)/%$(PRG_SUFFIX), $(PRG_SRC))

all : $(LOBJ) $(OBJ) $(PRG_OBJ) $(PRG)

%.o : %.c %.h
	$(CC) $(CFLAGS) $(INC) -o $@ -c $<

$(ODIR)/%.o : %.c %.h
	mkdir -p $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ -c $<

$(ODIR)/%.o : %.c
	mkdir -p $(ODIR)
	$(CC) $(CFLAGS) $(INC) -o $@ -c $<

$(PRG) : $(OBJ) $(PRG_OBJ)
	mkdir -p $(EDIR)
	$(CC) $(LDFLAGS) $(LOBJ) $^ -o $@

.PHONY : clean
clean :
	rm $(OBJ) $(PRG_OBJ) $(PRG)
