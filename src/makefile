#.PHONY: all clean

PROG=DenPro
CC=g++
#CC=icpc
COPT=-c -D_NO_HASHTABLE #-Wall
LOPT=-lz
SRC=$(PROG).cpp common.cpp Data.cpp TxtFile.cpp Distrib.cpp
HDR=$(PROG).h def.h common.h Data.h TxtFile.h Distrib.h
#SRC=*.cpp
#HDR=*.h
OBJ=$(SRC:.cpp=.o)
EXEC=$(PROG)

all: $(HDR) $(SRC) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(LOPT) $(OBJ) -o $@
	cp $@ ..

.cpp.o:
	$(CC) $(COPT) $< -o $@

#clean:
#	rm *o
