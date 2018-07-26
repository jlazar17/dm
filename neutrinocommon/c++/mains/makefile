# Compiler
GCC=g++

# Directories
LIBDIR        = ../lib
INCDIR        = ../inc
SUINCDIR      = ../inc/SU_inc
BINDIR        = ../bin

# other libraries locations
CUDAINC		 = /Developer/NVIDIA/CUDA-5.0/include
CUDALIB       = /Developer/NVIDIA/CUDA-5.0/lib
NUSIGMALIB    = ../../linked/nusigma/lib
GFORTRANLIB   = /usr/local/Cellar/gfortran/4.8.2/gfortran/lib 
#GFORTRANLIB   = /usr/local/gfortran/lib
GSLLIB   = /Users/carguelles/Programs/gsl/build/lib
NUSIGMA       = -lnusigma -lpythia -lxcm -lgfortran
CUDA          = -lcublas -lcudart

CFLAGS        = -Wall -g -O3
ifdef GPU_ENABLE
	INCFLAGS      = -I$(INCDIR) -I$(SUINCDIR) -I$(CUDAINC)
	LIBFLAGS      = -L$(LIBDIR) -L$(NUSIGMALIB) -L$(CUDALIB) -L$(GFORTRANLIB) -lneutrinocommon $(NUSIGMA) $(CUDA)
	LDFLAGS       = -lm -lgsl -lgslcblas -lstdc++ -lpthread
else
	INCFLAGS      = -I$(INCDIR) -I$(SUINCDIR) -I/Users/carguelles/Programs/gsl/build/include
	LIBFLAGS      = -L$(LIBDIR) -L$(NUSIGMALIB) -L$(GSLLIB) -L$(GFORTRANLIB) -lneutrinocommon $(NUSIGMA)
	LDFLAGS       = -lm -lgsl -lgslcblas -lstdc++ -lpthread
endif

all: demo test testinterface testsualg process process_int

demo: $(LIBDIR)/libneutrinocommon.a demo.cpp
	$(GCC) $(CFLAGS) $(INCFLAGS) -o $(BINDIR)/demo.exe demo.cpp $(LIBFLAGS) $(LDFLAGS)
	
test: $(LIBDIR)/libneutrinocommon.a test.cpp
	$(GCC) $(CFLAGS) $(INCFLAGS) -o $(BINDIR)/test.exe test.cpp $(LIBFLAGS) $(LDFLAGS)
	
process: $(LIBDIR)/libneutrinocommon.a process.cpp
	$(GCC) $(CFLAGS) $(INCFLAGS) -o $(BINDIR)/process.exe process.cpp $(LIBFLAGS) $(LDFLAGS)
	
process_int: $(LIBDIR)/libneutrinocommon.a process_int.cpp
	$(GCC) $(CFLAGS) $(INCFLAGS) -o $(BINDIR)/process_int.exe process_int.cpp $(LIBFLAGS) $(LDFLAGS)
	
testinterface: $(LIBDIR)/libneutrinocommon.a test_interface.cpp
	$(GCC) $(CFLAGS) $(INCFLAGS) -o $(BINDIR)/test_interface.exe test_interface.cpp $(LIBFLAGS) $(LDFLAGS)
	
testgputools: $(LIBDIR)/libneutrinocommon.a test_gputools.cpp
	$(GCC) $(CFLAGS) $(INCFLAGS) -o $(BINDIR)/test_gputools.exe test_gputools.cpp $(LIBFLAGS) $(LDFLAGS)
	
testsualg: $(LIBDIR)/libneutrinocommon.a test_sualg.cpp
	$(GCC) $(CFLAGS) $(INCFLAGS) -o $(BINDIR)/test_sualg.exe test_sualg.cpp $(LIBFLAGS) $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(BINDIR)/*
