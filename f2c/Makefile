
.PHONY : clean default

CFLAGS = -W -Wall -O3 -std=c99
LDLIBS = -lm

default : dependencies.pdf MyTestProgram TestVF

# Checked:
#   enorm()
#   lmdif1()
#   fdjac2()
#   qrfac()
#   qrsolve() -- but: recheck!

MyTestProgram : MyTestProgram.c lmdif1.c lmdif.c enorm.c qrfac.c lmpar.c qrsolv.c fdjac2.c minpack_c.h

TestVF : TestVF.c lmdif1.c lmdif.c enorm.c qrfac.c lmpar.c qrsolv.c fdjac2.c minpack_c.h

dependencies.pdf : dependencies.dot
	dot -Tpdf dependencies.dot -o $@
	gs -o $@.temp -sDEVICE=pdfwrite -sPAPERSIZE=a4 -dFIXEDMEDIA -dPDFFitPage -dCompatibilityLevel=1.4 $@
	mv $@.temp $@

clean :
	$(RM) *.o MyTestProgram TestVF dependencies.pdf
