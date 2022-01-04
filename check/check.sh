#!/bin/bash -ex
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#*                                                                           *#
#*                  This file is part of the class library                   *#
#*       SoPlex --- the Sequential object-oriented simPlex.                  *#
#*                                                                           *#
#*    Copyright (C) 1996-2022 Konrad-Zuse-Zentrum                            *#
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

SOPLEX_BOOST_SUPPORT="$(../${EXECUTABLE} --solvemode=2 |& grep 'without linking boost is not supported'|head -n 1)"
if [[ "${SOPLEX_BOOST_SUPPORT}" == "Using rational methods without linking boost is not supported" ]]
then
    SETTINGSLIST=(default devex steep)
else
    SETTINGSLIST=(default devex steep exact)
fi

if ! test -f ../settings/default.set
then
    touch ../settings/default.set
fi

# Solve with the different settings
for SETTINGS in ${SETTINGSLIST[@]}
do
    ./test.sh quick "${EXECUTABLE}" "${SETTINGS}" 60 "${OUTPUTDIR}"
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
