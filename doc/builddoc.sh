#!/bin/bash -e

### START SHELL TUTORIAL

DOXYFILE=${1:-soplex.dxy}

# build a fresh version of SoPlex
make -j clean -C ../
make -j -C ../

### PARAMETER FILE CREATION

../bin/soplex --saveset=parameters.set

cd inc
#parser.py now generated faq.inc as well
#python3 parser.py --linkext $HTML_FILE_EXTENSION  && php localfaq.php > faq.inc
./parser.py --linkext ${HTML_FILE_EXTENSION:=html}
cd ..

### FINISHED FAQ GENERATION

# finally build the soplex documentation
doxygen $DOXYFILE
