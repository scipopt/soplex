#!/usr/bin/env bash

# Solves a number of standard settings with a short testset

BINFILE=$1
BINNAME=`basename $BINFILE`

SETTINGSLIST=(default devex steep)

# Solve with the different settings
for SETTINGS in ${SETTINGSLIST[@]}
do
    ./test.sh quick $BINFILE $SETTINGS 60
done

echo
echo 'Results:'
for SETTINGS in ${SETTINGSLIST[@]}
do
    echo 'check/results/check.quick.'$BINNAME'.'$SETTINGS'.res'
done


# Evalulate the results
