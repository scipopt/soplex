#!/bin/bash -ex
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#*                                                                           *#
#*                  This file is part of the class library                   *#
#*       SoPlex --- the Sequential object-oriented simPlex.                  *#
#*                                                                           *#
#*    Copyright (C) 1996-2021 Konrad-Zuse-Zentrum                            *#
#*                            fuer Informationstechnik Berlin                *#
#*                                                                           *#
#*  SoPlex is distributed under the terms of the ZIB Academic Licence.       *#
#*                                                                           *#
#*  You should have received a copy of the ZIB Academic License              *#
#*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  *#
#*                                                                           *#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#

# Solves a number of standard settings with a short testset
# Called by 'make check' from soplex root

EXECUTABLE="${1}"
BINNAME=$(basename "${EXECUTABLE}")
HOST=$(uname -n | sed 's/\(.zib.de\)//g')
BINNAME="${BINNAME}.${HOST}"

OUTPUTDIR=results/quick

SETTINGSLIST=(default devex steep exact)

if ! test -f ../settings/default.set
then
    touch ../settings/default.set
fi

# Solve with the different settings
for SETTINGS in ${SETTINGSLIST[@]}
do
    ./test.sh quick "${EXECUTABLE}" "${SETTINGS}" 60 "${OUTPUTDIR}" 0
done

echo
echo 'Summary:'
for SETTINGS in ${SETTINGSLIST[@]}
do
    echo
    grep 'Results' -A1 ${OUTPUTDIR}'/check.quick.'${BINNAME}'.'${SETTINGS}'.res'
    echo 'check/'${OUTPUTDIR}'/check.quick.'${BINNAME}'.'${SETTINGS}'.res'
done

# Evalulate the results
