includedir=$(prefix)/include/Tauola
lib_LTLIBRARIES = libTauolaCxxInterface.la

libTauolaCxxInterface_la_SOURCES =            \
eventRecordInterfaces/TauolaHEPEVTEvent.cxx    \
eventRecordInterfaces/TauolaHEPEVTParticle.cxx \
tauolaCInterfaces/Tauola.cxx                  \
tauolaCInterfaces/TauolaParticle.cxx          \
tauolaCInterfaces/TauolaParticlePair.cxx      \
tauolaCInterfaces/TauolaEvent.cxx             \
tauolaCInterfaces/DecayList.cxx               \
tauolaFortranInterfaces/f_FilHep.cxx          \
tauolaFortranInterfaces/f_Init.cxx            \
tauolaFortranInterfaces/f_Decay.cxx           \
tauolaFortranInterfaces/tauola_extras.f       \
utilities/Log.cxx                             \
utilities/Plots.cxx

libTauolaCxxInterface_la_LIBADD = ../tauola-fortran/libTauolaFortran.la

include_HEADERS =                            \
eventRecordInterfaces/TauolaHEPEVTParticle.h  \
eventRecordInterfaces/TauolaHEPEVTEvent.h     \
tauolaCInterfaces/Tauola.h                   \
tauolaCInterfaces/TauolaParticle.h           \
tauolaCInterfaces/TauolaParticlePair.h       \
tauolaCInterfaces/TauolaEvent.h              \
tauolaCInterfaces/DecayList.h                \
tauolaFortranInterfaces/f_FilHep.h           \
tauolaFortranInterfaces/f_Init.h             \
tauolaFortranInterfaces/f_Decay.h            \
tauolaFortranInterfaces/f_Variables.h        \
utilities/Log.h                              \
utilities/Plots.h

if HAS_HEPMC
  libTauolaCxxInterface_la_SOURCES += \
  eventRecordInterfaces/TauolaHepMCEvent.cxx    \
  eventRecordInterfaces/TauolaHepMCParticle.cxx

  include_HEADERS += \
  eventRecordInterfaces/TauolaHepMCParticle.h  \
  eventRecordInterfaces/TauolaHepMCEvent.h
endif


INCLUDES = -I$(top_srcdir)/src/tauolaFortranInterfaces \
           -I$(top_srcdir)/src/tauolaCInterfaces \
           -I$(top_srcdir)/src/utilities \
           -I$(HEPMC_DIR)/include

libTauolaCxxInterface_la_FFLAGS = -g -O2 -fPIC -fno-automatic -fno-backslash -ffixed-line-length-132
