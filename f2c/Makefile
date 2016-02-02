
.PHONY : clean default

CFLAGS = -W -Wall -O3 -std=c99
LDLIBS = -lm

default : dependencies.pdf MyTestProgram

# Checked:
#   enorm()
#   lmdif1()
#   fdjac2()
#   qrfac() -- but: recheck!

MyTestProgram : MyTestProgram.c lmdif1.c lmdif.c lmder1.c lmder.c enorm.c qrfac.c lmpar.c qrsolv.c fdjac2_2nd.c minpack_c.h

dependencies.pdf : dependencies.dot
	dot -Tpdf dependencies.dot -o $@
	gs -o $@.temp -sDEVICE=pdfwrite -sPAPERSIZE=a4 -dFIXEDMEDIA -dPDFFitPage -dCompatibilityLevel=1.4 $@
	mv $@.temp $@

clean :
	$(RM) *.o MyTestProgram dependencies.pdf
