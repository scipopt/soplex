#!/bin/bash

### START SHELL TUTORIAL

if [ -n "$1" ]
then
    DOXYFILE=$1
else
    DOXYFILE=soplex.dxy
fi

# build a fresh version of SoPlex
make -j clean -C ../
make -j -C ../

### PARAMETER FILE CREATION

../bin/soplex --saveset=parameters.set

if [ -z "$HTML_FILE_EXTENSION" ]
then
    HTML_FILE_EXTENSION=html
fi

cd inc
python3 parser.py --linkext $HTML_FILE_EXTENSION  && php localfaq.php > faq.inc
cd ../

### FINISHED FAQ GENERATION

# finally build the soplex documentation
doxygen $DOXYFILE
