##############################################################################
# In order to compile and run TAUOLA without autoconfig                      #
# copy this file into  main TAUOLA directory and rename it  'make.inc'       #
# you will need to  adjust it to your platform                               #
# WARNING: 'make.inc' will be deleted if 'make Clean' is executed            #
#           use 'make clean'                                                 #
##############################################################################


###############################
# 1) paths                    #
###############################

# Absolute path to HepMC -- essential.
# Note: since configuration is not used, for HepMC 2.03 or older execute
#      './platform/to-HepMC-2.03.sh' script to convert the interface

HEPMCLOCATION =

# To compile examples (to run see below)

PYTHIALOCATION   = 
MCTESTERLOCATION = 

# Absolute path to installation directory (will be used if you 'make install')

PREFIX        = /usr/local/ 


###############################
# 2) compiler linker etc      #
###############################

# Fortran compiler and its flags. If you have g++ version 4.1 or later
# use gfortran, otherwise try g77 or f77. (check the g++ version by typing
# 'g++ --version'. Flags for fortran compiler normally should work unchanged

F77      = gfortran
F77FLAGS = -g -O2 -fPIC -fno-automatic -fno-backslash -ffixed-line-length-132
FFLAGS   = $(F77FLAGS)

# Linker. 'ld' works fine in most cases, but to ensure it's called properly,
# The fortran compiler is used to invoke the linker. These variables
# should normally be left intact

LD       = $(F77)
LDFLAGS  = -g
SOFLAGS  = -shared

RANLIB   = ranlib
AR       = ar
DIFF    = diff -ib

# Preprocessor. The best one is of course 'cpp', but if not accessible
# try 'g++ -E' or 'gcc -E'. All three flags are necessary

CPP     = cpp
CPPFLAGS= -C -P -traditional-cpp

# C++ compiler and its flags. If other than default version is needed,
# try 'g++-34', 'g++-41' or 'g++34', 'g++41' etc. (if they exist, which
# can be checked by typing the appropriate name in console)

CC      = g++
CFLAGS  = -g -O2 -fPIC -Wall -Wno-write-strings

# C compiler and its flags. Not used by TAUOLA, but needed by SANC library.
# Behavies in the same way as C++ compiler.

GCC     = gcc
GCCFLAGS= -O3 -fomit-frame-pointer -DHAVE_UNDERSCORE


###############################
# 3) for debugging            #
###############################

# Uncomment the line below, compile with debug mode turned on

#DEBUG = -D"_LOG_DEBUG_MODE_"


#################################
# 4) commands to run examples   #
#################################

##############################################################################
# invoke commands                                                            #
#    export PYTHIA8DATA=<path_to_pythia>/xmldoc                              #
#    export PATH=<path_to_root_binaries>:${PATH}                             #
#    export MCTESTERLOCATION=<path_to_MC_TESTER>                             #
#    touch  examples/make.inc                                                # 
##############################################################################


###############################
# 5) some info                #
###############################

##############################################################################
#                                                                            #
# PYTHIA 8.1: PYTHIALOCATION must be set to Pythia main directory. Pythia    #
#             must be compiled with HepMC. Global variable PYTHIA8DATA must  #
#             be set. Try:                                                   #
#             export PYTHIA8DATA=<path_to_pythia>/xmldoc                     #
#             <path_to_pythia> same as PYTHIALOCATION                        #
# MC-TESTER:  MCTESTERLOCATION must be set to MC-TESTER main directory.      #
#             It must be also set as global variable, for 'examples/testing' #
#             to work. Try:                                                  #
#             export MCTESTERLOCATION=<path_to_MC_TESTER>                    #
#             MC-TESTER must be compiled with HepMC                          #
#             ROOT must also be installed, as its required by MC-TESTER. Its #
#             binaries must be located in $PATH global variable. Check it by #
#             executing 'root' in console. If root is not set, try:          #
#             export PATH=<path_to_root_binaries>:${PATH}                    #
#             <path_to_MC_TESTER> same as MCTESTERLOCATION                   #
# WARNING:    Check that all is compiled/linked with same versions of        #
#             compilers, HepMC etc.                                          #
##############################################################################

