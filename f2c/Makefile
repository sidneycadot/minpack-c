
.PHONY : clean

CFLAGS = -W -Wall -O3 -std=c99
LDLIBS = -lm

C_FILES = $(MINPACK_FILES:%=%.c)


MyTestProgram : MyTestProgram.c lmdif1.c lmdif.c lmder1.c lmder.c enorm.c qrfac.c lmpar.c qrsolv.c fdjac2.c minpack_c.h dpmpar.h

clean :
	$(RM) *.o MyTestProgram
