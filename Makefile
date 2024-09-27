#----------------------------------------------------------------------------
#        makefile 
# type 'make' to compile the programs
#----------------------------------------------------------------------------

CC        = gcc

#----------------------------------------------------------------------------
# Project specific path defintions.
#----------------------------------------------------------------------------

PROJDIR    = .
SRC        = $(PROJDIR)/src
BIN        = $(PROJDIR)/bin

all : findframe importmats multiplymats pointmats makeassembly frac2orth \
	movecoords cif2pdb

findframe : $(SRC)/findframe.c 
	$(CC) -g -o $(BIN)/$@  $(SRC)/findframe.c  -lm

importmats : $(SRC)/importmats.c
	$(CC) -g -o $(BIN)/$@  $(SRC)/importmats.c -lm

multiplymats : $(SRC)/multiplymats.c
	$(CC) -g -o $(BIN)/$@  $(SRC)/multiplymats.c -lm

pointmats : $(SRC)/pointmats.c
	$(CC) -g -o $(BIN)/$@  $(SRC)/pointmats.c -lm

makeassembly : $(SRC)/makeassembly.c
	$(CC) -g -o $(BIN)/$@  $(SRC)/makeassembly.c -lm

frac2orth : $(SRC)/fracorth.c
	$(CC) -g -o $(BIN)/$@  $(SRC)/fracorth.c -lm

movecoords : $(SRC)/movecoords.c
	$(CC) -g -o $(BIN)/$@  $(SRC)/movecoords.c -lm

cif2pdb : $(SRC)/cif2pdb.c
	$(CC) -g -o $(BIN)/$@  $(SRC)/cif2pdb.c -lm

invertmat : $(SRC)/invertmat.c
        $(CC) -g -o $(BIN)/$@  $(SRC)/invertmat.c -lm

.PHONY: clean
clean :
	@cd $(BIN) ; rm -f findframe importmats multiplymats pointmats  \
	makeassembly frac2orth movecoords cif2pdb

