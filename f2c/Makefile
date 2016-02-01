
.PHONY : clean

CFLAGS = -W -Wall -O3 -std=c99
LDLIBS = -lm

C_FILES = $(MINPACK_FILES:%=%.c)


MyTestProgram : MyTestProgram.c lmdif1.c lmdif.c lmder.c enorm.c dpmpar.c qrfac.c lmpar.c qrsolv.c fdjac2.c

clean :
	$(RM) *.o MyTestProgram
