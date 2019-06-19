
CC := gcc-9
# CC := gcc
CFLAGS := -Wall -Wshadow -Winline -O3 -std=c11 -fopenmp

CXX := g++-9
# CXX := g++
CXXFLAGS := -Wall -Wshadow -Winline -O3 -std=c++11

INC := -I../../mtwister/ -I../../arralloc/
LDFLAGS := -lm -lgomp -larmadillo

PRG_SUFFIX_FLAG := 0

PRG_SUFFIX =

ODIR := ../bin/obj/
EDIR := ../bin/exe/

SRC_C := cell.c phase_field_model.c array.c image.c shape.c dump.c dump_cm.c dump_bulk_cm.c dump_field.c dump_cell_field.c dump_gyration.c dump_shape.c random_mtwister.c
OBJ_C := $(patsubst %.c, $(ODIR)/%.o, $(SRC_C))
SRC_CXX := filter.cpp sgolay_filter.cpp
OBJ_CXX := $(patsubst %.cpp, $(ODIR)/%.o, $(SRC_CXX))
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

$(ODIR)/%.o : %.cpp %.h
	mkdir -p $(ODIR)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ -c $<

$(ODIR)/%.o : %.cpp
	mkdir -p $(ODIR)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ -c $<

$(PRG) : $(OBJ_C) $(OBJ_CXX) $(PRG_OBJ)
	mkdir -p $(EDIR)
	$(CXX) $(LDFLAGS) $(LOBJ) $^ -o $@

.PHONY : clean
clean :
	rm -rf $(OBJ_C) $(OBJ_CXX) $(PRG_OBJ) $(PRG)
