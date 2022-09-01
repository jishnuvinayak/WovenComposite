##############################################################
# Makefile for compiling Yarn UMAT.
# Ref : Adapted from Monolithic FE2 package.
##############################################################
ABACALL=abq2020
MAINSRC = UMAT.f90
TRG=bin/$(patsubst %.f90,%-std.o,${MAINSRC})
SRC := $(wildcard src/*.f90) $(wildcard src/*.INC)

build: ${TRG}

#Note CPATH is read as include path by ifort
${TRG}: src/${MAINSRC} ${SRC}
	cd src;  ${ABACALL} make library=${MAINSRC}
	mv -f src/$(notdir $@)  bin/
	rm -f src/libstandardU.so
	

