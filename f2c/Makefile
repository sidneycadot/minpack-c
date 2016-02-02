
.PHONY : clean default

CFLAGS = -W -Wall -O3 -std=c99
LDLIBS = -lm

default : dependencies.pdf MyTestProgram

MyTestProgram : MyTestProgram.c lmdif1.c lmdif.c lmder1.c lmder.c enorm.c qrfac.c lmpar.c qrsolv.c fdjac2.c minpack_c.h

dependencies.pdf : dependencies.dot
	dot -Tpdf dependencies.dot -o dependencies.pdf

clean :
	$(RM) *.o MyTestProgram dependencies.pdf
