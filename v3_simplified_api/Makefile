
.PHONY : clean default

CFLAGS = -W -Wall -O3 -std=c99
LDLIBS = -lm

default : TestVF TestVF_Performance

TestVF : TestVF.o lmdif1.o lmdif.o enorm.o qrfac.o lmpar.o qrsolv.o fdjac2.o

TestVF_Performance : TestVF_Performance.o VF_Table.o lmdif1.o lmdif.o enorm.o qrfac.o lmpar.o qrsolv.o fdjac2.o

TestVF.o : TestVF.c minpack_c.h

TestVF_Performance.o : TestVF.c minpack_c.h VF_Table.c

VF_Table.o : VF_Table.c VF_Table.h

lmdif1.o : lmdif1.c minpack_c.h
lmdif.o  : lmdif.c  minpack_c.h
enorm.o  : enorm.c  minpack_c.h
qrfac.o  : qrfac.c  minpack_c.h
lmpar.o  : lmpar.c  minpack_c.h
qrsolv.o : qrsolv.c minpack_c.h
fdjac2.o : fdjac2.c minpack_c.h

clean :
	$(RM) *.o TestVF TestVF_Performance
