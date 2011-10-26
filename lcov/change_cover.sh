#!/bin/sh

TEST_SUITE=coverage

# Determine base directory
BASE_DIR=`cd ..; pwd`
OBJ_DIR=obj/O.linux.x86.gnu.gcov.static/lib

export LD_LIBRARY_PATH=/opt/gcc-400/lib
export PATH=/opt/gcc-400/bin:$PATH

cd ..
make OPT=gcov change_exerciser

cd $OBJ_DIR
rm -f *.gcda
cd ../../..

cd check
../bin/exercise_LP_changes.linux.x86.gnu.gcov.static $TEST_SUITE.test
cd ..

cd lcov
lcov --capture -directory $BASE_DIR/$OBJ_DIR --output-file change_exerciser.info --test-name change_exerciser 
genhtml change_exerciser.info --output-directory change_exerciser --title "change_exerciser" --show-details