
.PHONY : clean default

CXXFLAGS = -W -Wall -O3

default : dependencies.pdf show_dpmpar

dependencies.pdf : dependencies.dot
	dot -Tpdf dependencies.dot -o $@
	#gs -o $@.temp -sDEVICE=pdfwrite -sPAPERSIZE=a4 -dFIXEDMEDIA -dPDFFitPage -dCompatibilityLevel=1.4 $@
	#mv $@.temp $@

show_dpmpar : show_dpmpar.cc

clean :
	$(RM) dependencies.pdf show_dpmpar
