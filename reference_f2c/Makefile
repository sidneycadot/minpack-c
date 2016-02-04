
.PHONY : clean default

CFLAGS = -W -Wall -O3 -std=c99
LDLIBS = -lm

default : TestVF

TestVF : TestVF.c lmdif1.c lmdif.c enorm.c qrfac.c lmpar.c qrsolv.c fdjac2.c dpmpar.c my_include.h

clean :
	$(RM) *.o TestVF
