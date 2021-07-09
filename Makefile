
CC=gcc
CFLAGS=-I.

DEPS = lattice.h plot.h gauss.h jacob.h over_relaxation.h charge_cage.h

OBJ = relaxation_method.o lattice.o plot.o gauss.o jacob.o over_relaxation.o charge_cage.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

relaxation_method: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)


.PHONY: clean

clean:
	rm -f *.o
