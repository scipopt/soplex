#!/bin/bash

### START SHELL TUTORIAL

# build a fresh version of SoPlex
make -j clean -C ../
make -j -C ../

### PARAMETER FILE CREATION

../bin/soplex --saveset=parameters.set

# finally build the soplex documentation
doxygen soplex.dxy
