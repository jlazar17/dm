exampledir = $(top_srcdir)/examples

### Standalone example ###
example_PROGRAMS = taumain_hepevt_example.exe

taumain_hepevt_example_exe_SOURCES  = taumain_hepevt_example.cxx

INCLUDES = -I$(prefix)/include

AM_LDFLAGS = -R $(prefix)/lib

LDADD =                          \
$(FLIBS)                         \
$(prefix)/lib/libTauolaCxxInterface.so  \
$(prefix)/lib/libTauolaFortran.so

### All other examples require HepMC ###
if HAS_HEPMC

AM_LDFLAGS += -R $(HEPMC_DIR)/lib

INCLUDES += -I$(HEPMC_DIR)/include
LDADD    += -L$(HEPMC_DIR)/lib -lHepMC

example_PROGRAMS += taumain_stand_alone_example.exe

taumain_stand_alone_example_exe_SOURCES  = taumain_stand_alone_example.cxx

### Taugun & taummk examples (require Pythia8) ###
if HAS_PYTHIA8

AM_LDFLAGS  += -R $(PYTHIA8_DIR)/lib/archive

INCLUDES += -I$(PYTHIA8_DIR)/include
LDADD    += -L$(PYTHIA8_DIR)/lib/archive -lpythia8 -llhapdfdummy -lhepmcinterface

example_PROGRAMS += single_tau_gun_example.exe taummk_pythia_example.exe

single_tau_gun_example_exe_SOURCES  = single_tau_gun_example.cxx
taummk_pythia_example_exe_SOURCES   = taummk_pythia_example.cxx

### Main example (requires Pythia8 and MC-Tester) ###
if HAS_MCTESTER

AM_LDFLAGS  += -R $(MCTESTER_DIR)/lib

INCLUDES += -I$(PYTHIA8_DIR)/include -I$(MCTESTER_DIR)/include -I$(ROOTINC)
LDADD    += -L$(MCTESTER_DIR)/lib -lHEPEvent -lHepMCEvent -lMCTester $(ROOTLIBS)

example_PROGRAMS += taumain_pythia_example.exe

taumain_pythia_example_exe_SOURCES  = taumain_pythia_example.cxx

endif # MC-TESTER

endif # PYTHIA8

endif # HepMC
