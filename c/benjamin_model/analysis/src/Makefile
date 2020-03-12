
CC := gcc
CFLAGS := -O3 -Wall -Wshadow -std=c11 -fopenmp

CXX := g++
CXXFLAGS := -O3 -Wall -Wshadow -std=c++11 -fopenmp
INC = ./lib/
LDFLAGS := -larmadillo -lgomp

PRG_SUFFIX_FLAG := 0

LDIR := ./lib
LODIR = ../bin/lib
ODIR := ../bin/obj
EDIR := ../bin/exe

LSRCC := $(wildcard $(LDIR)/*.c)
LOBJC := $(patsubst %.c, $(LODIR)/%.o, $(notdir $(LSRCC)))
LSRCCXX := $(wildcard $(LDIR)/*.cpp)
LOBJCXX := $(patsubst %.cpp, $(LODIR)/%.o, $(notdir $(LSRCCXX)))
SRC := $(wildcard *.cpp)
OBJ := $(patsubst %.cpp, $(ODIR)/%.o, $(SRC))
PRG := $(patsubst %.cpp,$(EDIR)/%,$(SRC))
PRG_SUFFIX=
BIN := $(patsubst %,%$(PRG_SUFFIX),$(PRG))


all : $(LOBJC) $(LOBJCXX) $(OBJ) $(BIN)

% : %.cpp
	make $(EDIR)/$@

$(ODIR)/%.o : %.cpp 
	mkdir -p $(ODIR)
	$(CXX) $(CXXFLAGS) -I$(INC) -o $@ -c $<

$(EDIR)/% : $(ODIR)/%.o $(LOBJC) $(LOBJCXX)
	mkdir -p $(EDIR)
	$(CXX) $(LDFLAGS) $^ -o $@

$(LODIR)/%.o : $(LDIR)/%.c
	mkdir -p $(LODIR)
	$(CC) $(CFLAGS) -o $@ -c $<

$(LODIR)/%.o : $(LDIR)/%.cpp
	mkdir -p $(LODIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.PHONY : clean
clean : 
	rm -rf $(EDIR) $(ODIR) $(LODIR)
