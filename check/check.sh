#!/usr/bin/env bash

# Solves a number of standard settings with a short testset

BINFILE=$1
BINNAME=`basename $BINFILE`
HOST=`uname -n | sed 's/\(.zib.de\)//g'`
BINNAME=$BINNAME.$HOST

RESDIR=results/quick

SETTINGSLIST=(default devex steep exact)

if ! test -f ../settings/default.set
then
    touch ../settings/default.set
fi

# Solve with the different settings
for SETTINGS in ${SETTINGSLIST[@]}
do
    ./test.sh quick $BINFILE $SETTINGS 60 $RESDIR
done

echo
echo 'Summary:'
for SETTINGS in ${SETTINGSLIST[@]}
do
    echo
    grep 'Results' -A1 $RESDIR'/check.quick.'$BINNAME'.'$SETTINGS'.res'
    echo 'check/'$RESDIR'/check.quick.'$BINNAME'.'$SETTINGS'.res'
done


# Evalulate the results
