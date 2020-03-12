
CC := gcc
# CC := gcc-9
CFLAGS := -Wall -Wshadow -Winline -O3 -std=c11 -fopenmp
#CFLAGS := -Wall -Wshadow -Winline -O3 -std=c11

CXX := g++
# CXX := g++-9
CXXFLAGS := -Wall -Wshadow -Winline -O3 -std=c++11

INC := -I../../mtwister/ -I../../arralloc/
LFLAGS := -lm -lgomp
#LFLAGS := -lm

PRG_SUFFIX_FLAG := 0
PRG_SUFFIX =

# 0 = no armadillo and 1 = with armadillo
HAS_ARMA := 1

#ODIR := ../bin/obj/
#EDIR := ../bin/exe/
ODIR := ../new_bin/obj/
EDIR := ../new_bin/exe/

SRC_C := cell.c phase_field_model.c array.c image.c overlap.c dump.c dump_cm.c dump_bulk_cm.c dump_field.c dump_cell_field.c dump_gyration.c dump_neighbour.c dump_energy.c dump_overlap.c dump_overlap_field.c dump_index_field.c random_mtwister.c
SRC_CXX := neighbour.cpp
ifeq ($(HAS_ARMA),1)
	SRC_C := $(SRC_C) shape.c dump_shape.c
	SRC_CXX := $(SRC_CXX) sgolay_filter.cpp
	LFLAGS := $(LFLAGS) -larmadillo
endif
OBJ_C := $(patsubst %.c, $(ODIR)/%.o, $(SRC_C))
OBJ_CXX := $(patsubst %.cpp, $(ODIR)/%.o, $(SRC_CXX))
PRG_SRC := run_phase_field_model.c
PRG_OBJ := $(patsubst %.c, $(ODIR)/%.o, $(PRG_SRC))
LSRC := ../../mtwister/mtwister.c ../../arralloc/arralloc.c
LOBJ := $(patsubst %.c, %.o, $(LSRC)) #../mtwister/mtwister.o ../arralloc/arralloc.o
PRG := $(patsubst %.c, $(EDIR)/%$(PRG_SUFFIX), $(PRG_SRC))

all : check_arma make_bin $(LOBJ) $(OBJ) $(PRG_OBJ) $(PRG)

%.o : %.c %.h
	$(CC) $(CFLAGS) $(INC) -o $@ -c $<

$(ODIR)/%.o : %.c %.h
	$(CC) $(CFLAGS) $(INC) -o $@ -c $<

$(ODIR)/%.o : %.c
	$(CC) $(CFLAGS) $(INC) -o $@ -c $<

$(ODIR)/%.o : %.cpp %.hpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ -c $<

$(ODIR)/%.o : %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INC) -o $@ -c $<

$(ODIR)/%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ -c $<

$(PRG) : $(OBJ_C) $(OBJ_CXX) $(PRG_OBJ)
	$(CXX) $(LFLAGS) $(LOBJ) $^ -o $@

.PHONY : clean check_arma make_bin

clean :
	rm -rf $(OBJ_C) $(OBJ_CXX) $(PRG_OBJ) $(PRG)

check_arma :
ifeq ($(HAS_ARMA),1)
	perl -pi -e 's/PF_HAS_ARMA [0,1]/PF_HAS_ARMA 1/g' dump.h
else
	perl -pi -e 's/PF_HAS_ARMA [0,1]/PF_HAS_ARMA 0/g' dump.h
endif

make_bin :
	mkdir -p $(ODIR)
	mkdir -p $(EDIR)
