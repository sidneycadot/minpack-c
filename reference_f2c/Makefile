
.PHONY : clean default

CFLAGS = -W -Wall -O3 -std=c99
LDLIBS = -lm

default : TestVF

TestVF : TestVF.o lmdif1.o lmdif.o enorm.o qrfac.o lmpar.o qrsolv.o fdjac2.o dpmpar.o

TestVF.o : TestVF.c my_include.h

lmdif1.o : lmdif1.c my_include.h
lmdif.o  : lmdif.c  my_include.h
enorm.o  : enorm.c  my_include.h
qrfac.o  : qrfac.c  my_include.h
lmpar.o  : lmpar.c  my_include.h
qrsolv.o : qrsolv.c my_include.h
fdjac2.o : fdjac2.c my_include.h
dpmpar.o : dpmpar.c my_include.h


clean :
	$(RM) *.o TestVF
