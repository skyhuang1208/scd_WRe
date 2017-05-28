# Specify which compiler to use:
CC= gcc

# Use CFLAGS to specify additional compiler options to use:
CFLAGS= -Wall -g -DHASH_DEBUG=1
#CFLAGS= -Wall -O3

# Define variables:
EXE_FILE= scdexe
OBJ_FILES= main.o hash.o rates.o rvgs.o RateMatrix.o
HEADER_FILES= uthash.h types.h constants.h rvgs.h RateMatrix.h

# Build executable:
exe : $(EXE_FILE)
$(EXE_FILE) : $(OBJ_FILES)
	$(CC) -o $(EXE_FILE) $(OBJ_FILES) -g -lm
	chmod ugo+x $(EXE_FILE)	# turn on execute permissions

# Rules:
.c.o:	$(HEADER_FILES)
	$(CC) $(CFLAGS) -c $*.c

clean:
	rm $(EXE_FILE) $(OBJ_FILES)

