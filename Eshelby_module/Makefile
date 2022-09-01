##############################################################
# Makefile for creating the Eshelby module wrapper.
##############################################################
ABACALL=abq2020
SRC = src/Eshelby_tensor.f90 src/Gauss_quadrature.f90 src/Eshelby_utilities.f90
TRG=Eshelby

${TRG}: ${SRC}
	${ABACALL} python -m numpy.f2py -m $@ -c $< 
	mv -f $@.so  bin/
	rm -f *.mod
