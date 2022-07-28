NHDIR = $(dir $(CURDIR))

# This needs to be set in order to find the febio .h files.
#FEBINC = $(HOME)/febio-2.0new/sdk/include
FEBINC = $(HOME)/codes/FEBio
# This needs to be set in order to find the febio library files.
#FEBLIB = $(HOME)/febio-2.0new/sdk/lib
FEBLIB = $(FEBINC)/cbuild/lib

lnx64: PLAT = lnx64
lnx32: PLAT = lnx32
osx: PLAT = osx
gcc: PLAT = gcc
lnx64clean: PLAT = lnx64
lnx32clean: PLAT = lnx32
osxclean: PLAT = osx
gccclean: PLAT = gcc

export PLAT
export FEBLIB
export NHDIR

# Assumes the Intel C++ compiler.  Use g++ for the gnu compiler and run "make gcc".
export CC = c++
export FLG = -O3 -fPIC -std=gnu++11

# Can be removed if not using the Intel compiler.
export LNKFLG = -static-intel
export INC = -I$(FEBINC)

lnx64 lnx32 gcc:
	( cd $(PLAT); $(MAKE) -f ../neohookeanpi.mk SFX="so" )

osx:
	( cd $(PLAT); $(MAKE) -f ../neohookeanpi.mk SFX="dylib" )

lnx64clean lnx32clean gccclean:
	( cd $(PLAT); $(MAKE) -f ../neohookeanpi.mk SFX="so" clean )

osxclean:
	( cd $(PLAT); $(MAKE) -f ../neohookeanpi.mk SFX="dylib" clean )

.PHONY: lnx64 lnx32 osx gcc
