
.PHONY : clean default

CFLAGS = -W -Wall -O3 -std=c99
LDLIBS = -lm

default : dependencies.pdf Test1 TestVF

# Checked:
#   enorm()
#   lmdif1()
#   fdjac2()
#   qrfac()
#   qrsolve() -- but: recheck!

Test1  : Test1.o  lmdif1.o lmdif.o enorm.o qrfac.o lmpar.o qrsolv.o fdjac2.o
TestVF : TestVF.o lmdif1.o lmdif.o enorm.o qrfac.o lmpar.o qrsolv.o fdjac2.o

Test1.o  : Test1.c  minpack_c.h
TestVF.o : TestVF.c minpack_c.h

lmdif1.o : lmdif1.c minpack_c.h
lmdif.o  : lmdif.c  minpack_c.h
enorm.o  : enorm.c  minpack_c.h
qrfac.o  : qrfac.c  minpack_c.h
lmpar.o  : lmpar.c  minpack_c.h
qrsolv.o : qrsolv.c minpack_c.h
fdjac2.o : fdjac2.c minpack_c.h
dpmpar.o : dpmpar.c minpack_c.h

dependencies.pdf : dependencies.dot
	dot -Tpdf dependencies.dot -o $@
	#gs -o $@.temp -sDEVICE=pdfwrite -sPAPERSIZE=a4 -dFIXEDMEDIA -dPDFFitPage -dCompatibilityLevel=1.4 $@
	#mv $@.temp $@

clean :
	$(RM) *.o Test1 TestVF dependencies.pdf
