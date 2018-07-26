PREFIX        = /tmp/TAUOLA.1.0.7/TAUOLA
HEPMCLOCATION = 

F77      = gfortran
F77FLAGS = -g -O2 -fPIC -fno-automatic -fno-backslash -ffixed-line-length-132
FFLAGS   = $(F77FLAGS)

LD       = gfortran
LDFLAGS  = -g
SOFLAGS  = -shared

RANLIB   = ranlib
AR       = ar
DIFF    = diff -ib

CPP     = cpp
CPPFLAGS= -C -P -traditional-cpp

CC      = g++
CFLAGS  = -g -O2 -fPIC -Wall -Wno-write-strings

GCC     = gcc
GCCFLAGS= -O3 -fomit-frame-pointer -DHAVE_UNDERSCORE

#Uncomment to compile with debug mode turned on
#DEBUG = -D"_LOG_DEBUG_MODE_"
