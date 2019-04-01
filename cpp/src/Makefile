
CXX := g++-8
CXXFLAGS := -O3 -Wall -Wshadow -std=c++11 -fopenmp
INC = -I./lib/ 
LDFLAGS := -lm -lomp

PRG_SUFFIX_FLAG := 0

LDIR := ./lib
LODIR = ../bin/lib
ODIR := ../bin/obj
EDIR := ../bin/exe

LSRC := $(wildcard $(LDIR)/*.cpp)
LOBJ := $(patsubst %.cpp, $(LODIR)/%.o, $(notdir $(LSRC)))
SRC := $(wildcard *.cpp)
OBJ := $(patsubst %.cpp, $(ODIR)/%.o, $(SRC))
PRG := $(EDIR)/RunPhaseFieldModel #$(patsubst %.cpp,$(EDIR)/%,$(SRC))
PRG_SUFFIX=
BIN := $(patsubst %,%$(PRG_SUFFIX),$(PRG))


all : $(LOBJ) $(OBJ) $(BIN)

% : %.cpp
	make $(EDIR)/$@

$(ODIR)/%.o : %.cpp 
	mkdir -p $(ODIR)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ -c $<

$(EDIR)/% : $(ODIR)/%.o $(LOBJ) $(OBJ)
	mkdir -p $(EDIR)
	$(CXX) $(LDFLAGS) $^ -o $@

$(LODIR)/%.o : $(LDIR)/%.cpp
	mkdir -p $(LODIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.PHONY : clean
clean : 
	rm -rf $(EDIR) $(ODIR) $(LODIR)
