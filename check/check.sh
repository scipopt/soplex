#!/usr/bin/env bash

# Solves a number of standard settings with a short testset

BINFILE=$1
BINNAME=`basename $BINFILE`

RESDIR=results/quick

SETTINGSLIST=(default devex steep)

# Solve with the different settings
for SETTINGS in ${SETTINGSLIST[@]}
do
    ./test.sh quick $BINFILE $SETTINGS 60 $RESDIR
done

echo
echo 'Summary:'
echo
for SETTINGS in ${SETTINGSLIST[@]}
do
    echo 'check/'$RESDIR'/check.quick.'$BINNAME'.'$SETTINGS'.res'
    grep 'Results:' -A1 $RESDIR'/check.quick.'$BINNAME'.'$SETTINGS'.res'
    echo
done


# Evalulate the results
