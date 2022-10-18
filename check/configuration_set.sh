#!/bin/bash
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

# configures environment variables for test runs both on the cluster and locally
# to be invoked inside a check(_cluster.sh) script
# This script cancels the process if required variables are not correctly set

# new environment variables defined by this script:
#    SOPLEXPATH - absolute path to invocation working directory
#    EXECUTABLE - absolute path to executable
#    FULLTSTNAME - path to .test file
#    SOLUFILE - .solu file for this test set, for parsing optimal solution values
#    SETTINGSFILE - absolute path to settings file

# get current SOPLEX path
SOPLEXPATH=$(pwd -P)

EXECUTABLE="${SOPLEXPATH}/../${EXECUTABLE}"

# search for test file in check/instancedata/testsets and in check/testset
if test -f "instancedata/testsets/${TSTNAME}.test"
then
    FULLTSTNAME="${SOPLEXPATH}/instancedata/testsets/${TSTNAME}.test"
elif test -f "testset/${TSTNAME}.test"
then
    FULLTSTNAME="${SOPLEXPATH}/testset/${TSTNAME}.test"
else
    echo "Skipping test: no ${TSTNAME}.test file found in testset/ or instancedata/testsets/"
    exit 1
fi

SETTINGSFILE="${SOPLEXPATH}/../settings/${SETTINGS}.set"
# Abort if files are missing
if ! test -f "${SETTINGSFILE}"
then
    if [ "${SETTINGSFILE}" == "${SOPLEXPATH}/../settings/default.set" ];
    then
        touch "${SETTINGSFILE}"
    else
        echo "Settings file not found: ${SETTINGSFILE}"
        exit 1
    fi
fi

# look for solufiles under the name of the test, the name of the test with everything after the first "_" or "-" stripped, and all;
# prefer more specific solufile names over general ones and the instance database solufiles over those in testset/
SOLUFILE=""
for f in "${TSTNAME}" ${TSTNAME%%_*} ${TSTNAME%%-*} all
do
    for d in instancedata/testsets testset
    do
        if test -f "${SOPLEXPATH}/${d}/${f}.solu"
        then
            SOLUFILE="${SOPLEXPATH}/${d}/${f}.solu"
            break
        fi
    done
    if ! test "${SOLUFILE}" = ""
    then
        break
    fi
done

