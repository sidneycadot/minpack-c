
.PHONY : clean

CFLAGS = -W -Wall -O3
LDFLAGS =  -lm

C_FILES = $(MINPACK_FILES:%=%.c)


MyTestProgram : MyTestProgram.c lmdif1.c lmder.c enorm.c dpmpar.c qrfac.c lmpar.c qrsolv.c

clean :
	$(RM) *.o MyTestProgram
