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
echo 'Results:'
for SETTINGS in ${SETTINGSLIST[@]}
do
    echo 'check/'$RESDIR'/check.quick.'$BINNAME'.'$SETTINGS'.res'
done


# Evalulate the results
