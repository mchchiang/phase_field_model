
CC := gcc # gcc-8
CFLAGS := -Wall -Wshadow -O3 -std=c11 -fopenmp
INC := -I../mtwister/ -I../arralloc/
LDFLAGS := -lm -lgomp

PRG_SUFFIX_FLAG := 0

PRG_SUFFIX =

SRC := cell.c phase_field_model.c array.c
OBJ := $(patsubst %.c, %.o, $(SRC))
PRG_SRC := run_phase_field_model.c
PRG_OBJ := $(patsubst %.c, %.o, $(PRG_SRC))
LSRC := ../mtwister/mtwister.c ../arralloc/arralloc.c
LOBJ := $(patsubst %.c, %.o, $(LSRC)) #../mtwister/mtwister.o ../arralloc/arralloc.o
PRG := $(patsubst %.c, %$(PRG_SUFFIX), $(PRG_SRC))

all : $(LOBJ) $(OBJ) $(PRG_OBJ) $(PRG)

%.o : %.c
	$(CC) $(CFLAGS) $(INC) -o $@ -c $<

$(PRG) : $(OBJ) $(PRG_OBJ)
	$(CC) $(LDFLAGS) $(LOBJ) $^ -o $@


.PHONY : clean
clean :
	rm $(OBJ) $(PRG_OBJ) $(PRG)
