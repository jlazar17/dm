# Compiler
GCC=g++
NVCC=nvcc
FC=gfortran

PROOT=/Users/carguelles/Programs/neutrinocommons/neutrinocommon/c++/

# Directories
LIBDIR        = ../lib
INCDIR        = ../inc
SUINCDIR      = ../inc/SU_inc
NUSIGMAINC    = ../../linked/nusigma/inc
NUSIGMALIB    = ../../linked/nusigma/lib
# change to your cuda directories
CUDAINC	      = /Developer/NVIDIA/CUDA-5.0/include
CUDALIB       = /Developer/NVIDIA/CUDA-5.0/lib
# cuda flags
CUFLAGS		  =  -lcublas -lcudart
INCCUFLAGS    = -I$(INCDIR)
# C flags
ifdef GPU_ENABLE
	CFLAGS        = -Wall -O3 -fPIC $(CUFLAGS)
	INCCFLAGS     = -I$(INCDIR) -I$(SUINCDIR) -I$(CUDAINC)
else
	CFLAGS        = -Wall -O3 -fPIC
	INCCFLAGS     = -I$(INCDIR) -I$(SUINCDIR) -I/Users/carguelles/Programs/gsl/build/include
endif

LIBFLAGS      = 
LDFLAGS       = -lm -lgsl -lgslcblas -lstdc++ -lpthread -lnusigma
LDFLAGS       += -L$(NUSIGMALIB) -L/Users/carguelles/Programs/gsl/build/lib

# F flags
FFLAGS		  = -O3 -fPIC
INCFFLAGS	  = -I$(NUSIGMAINC) -L$(NUSIGMALIB)


# Project files
LIBNAME = libneutrinocommon.a
FOBJ	= nusigma_c++.fo

ifdef GPU_ENABLE
	COBJ    = body.o neuosc.o physconst.o tools.o xsections.o neurho.o interface.o taudecay.o gpu_tools.o neualg.o SUNalg.o
	HEADER  = $(INCDIR)/body.h $(INCDIR)/neuosc.h $(INCDIR)/physconst.h $(INCDIR)/global.h $(INCDIR)/tools.h $(INCDIR)/xsections.h $(INCDIR)/neurho.h $(INCDIR)/interface.h $(INCDIR)/taudecay.h $(INCDIR)/gpu_tools.h $(INCDIR)/neualg.h $(INCDIR)/SUNalg.h
else
	COBJ    = body.o neuosc.o physconst.o tools.o xsections.o neurho.o interface.o taudecay.o neualg.o SUNalg.o
	HEADER  = $(INCDIR)/body.h $(INCDIR)/neuosc.h $(INCDIR)/physconst.h $(INCDIR)/global.h $(INCDIR)/tools.h $(INCDIR)/xsections.h $(INCDIR)/neurho.h $(INCDIR)/interface.h $(INCDIR)/taudecay.h $(INCDIR)/neualg.h $(INCDIR)/SUNalg.h
endif

# Compilation rules
lib libneutrinocommon $(LIBNAME): $(COBJ) $(FOBJ) makefile $(HEADER)
	ar -crv $(LIBDIR)/$(LIBNAME) $(COBJ) $(FOBJ)
	ranlib $(LIBDIR)/$(LIBNAME)
	$(GCC) -dynamiclib -install_name $(PROOT)/lib/libneutrinocommon.dylib $(LDFLAGS) $(COBJ) $(FOBJ) -o $(LIBDIR)/libneutrinocommon.dylib 

%.fo : %.f 
	$(FC) -c $(FFLAGS) $(INCFFLAGS) -o $(FOBJ) $<

%.o : %.cpp $(HEADER) $(FOBJ)
	$(GCC) -c $(CFLAGS) $(INCCFLAGS) $<

.PHONY: clean
clean:
	rm -f $(LIBDIR)/$(LIBNAME) $(COBJ) $(FOBJ) 
	rm -f *.o *.fo
