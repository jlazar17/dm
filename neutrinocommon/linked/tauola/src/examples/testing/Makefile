include ../make.inc                                                   
include $(TAUOLALOCATION)/make.inc                                       

ROOTLIB       :=  $(shell root-config --glibs)
ROOTINC       :=  $(shell root-config --cflags)
ROOTLIBPATH   :=  $(shell root-config --libdir)
PYTHIA_LIBS    = -L$(PYTHIALOCATION)/lib/archive -lpythia8 -llhapdfdummy -lhepmcinterface
MC_TESTER_LIBS = -L$(MCTESTERLOCATION)/lib -lHEPEvent -lHepMCEvent -lMCTester
HEPMC_LIBS     = -L$(HEPMCLOCATION)/lib -lHepMC
TAUOLA_LIBS    = -L$(TAUOLALOCATION)/lib -lTauolaCxxInterface -lTauolaFortran

all:
ifneq ($(HEPMCLOCATION), )
	@make EW-PLOTS/plots.exe
ifneq ($(PYTHIALOCATION), )
ifneq ($(MCTESTERLOCATION), )
	@make tauola_test.exe
endif
endif
endif
	@echo "##################################################################"
	@echo " Available programs:"
	@if test -e EW-PLOTS/plots.exe; then \
	  echo " ./EW-PLOTS/plots.exe"; \
	  echo "   to execute follow instructions in EW-PLOTS/README"; \
	fi
	@if test -e tauola_test.exe; then \
	  echo " ./tauola_test.exe"; \
	  echo "   to execute all tauola_test-s: 'make run',"; \
	  echo "   for single test: cd to chosen sub-directory and 'make' "; else \
	  echo ""; \
	  echo " WARNING: Main program used for tests - 'tauola_test.exe' "; \
          echo "          requires HepMC, Pythia8 and MC-TESTER."; \
	fi
	@echo "##################################################################"

%.o: %.c
	$(CC) $(CFLAGS) -I$(TAUOLALOCATION)/include -I$(PYTHIALOCATION)/include -I$(MCTESTERLOCATION)/include \
	-I$(HEPMCLOCATION)/include $(ROOTINC) -c $< -o $@

%.exe: %.o
	$(CC) -o $@ $(LDFLAGS) $< $(ROOTLIB) \
	$(PYTHIA_LIBS) $(MC_TESTER_LIBS) $(HEPMC_LIBS) $(TAUOLA_LIBS) \
	-Wl,-rpath,$(TAUOLALOCATION)/lib \
	-Wl,-rpath,$(ROOTLIBPATH) \
	-Wl,-rpath,$(HEPMCLOCATION)/lib \
	-Wl,-rpath,$(PYTHIALOCATION)/lib \
	-Wl,-rpath,$(MCTESTERLOCATION)/lib

# EW-PLOTS/plots.c - only HepMC required
EW-PLOTS/plots.exe: EW-PLOTS/plots.c
	$(CC) $(CFLAGS) -I$(TAUOLALOCATION)/include -I$(HEPMCLOCATION)/include -c $< -o $(<:.c=.o)
	$(CC) -o $@ $(LDFLAGS) $(<:.c=.o) $(HEPMC_LIBS) $(TAUOLA_LIBS) \
	-Wl,-rpath,$(TAUOLALOCATION)/lib \
	-Wl,-rpath,$(HEPMCLOCATION)/lib

make.inc:
	@echo ""
	@echo "Please execute ./configure in main examples directory first!"
	@echo ""
	@false

$(TAUOLALOCATION)/make.inc:
	@echo ""
	@echo "Please compile the interface first!"
	@echo ""
	@false

run:
	make -C H
	make -C Z
	make -C W
	make -C tau
	make -C H+-


clean:
	rm -f *.o *~ *.exe
	rm -f EW-PLOTS/*.o EW-PLOTS/*~ EW-PLOTS/*.exe

clobber: clean
	make -C H clobber
	make -C Z clobber
	make -C W clobber
	make -C H+- clobber
	make -C H clobber
	make -C tau clobber
