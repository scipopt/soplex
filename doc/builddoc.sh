#!/bin/bash

### START SHELL TUTORIAL

DOXYFILE=soplex.dxy
if [ -z $1 ]; then DOXYFILE=$1; fi

# build a fresh version of SoPlex
make -j clean -C ../
make -j -C ../

### PARAMETER FILE CREATION

../bin/soplex --saveset=parameters.set

# finally build the soplex documentation
doxygen $DOXYFILE
