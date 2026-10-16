#!/bin/bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#*                                                                           *#
#*                  This file is part of the class library                   *#
#*       SoPlex --- the Sequential object-oriented simPlex.                  *#
#*                                                                           *#
#*  Copyright (c) 1996-2026 Zuse Institute Berlin (ZIB)                      *#
#*                                                                           *#
#*  Licensed under the Apache License, Version 2.0 (the "License");          *#
#*  you may not use this file except in compliance with the License.         *#
#*  You may obtain a copy of the License at                                  *#
#*                                                                           *#
#*      http://www.apache.org/licenses/LICENSE-2.0                           *#
#*                                                                           *#
#*  Unless required by applicable law or agreed to in writing, software      *#
#*  distributed under the License is distributed on an "AS IS" BASIS,        *#
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *#
#*  See the License for the specific language governing permissions and      *#
#*  limitations under the License.                                           *#
#*                                                                           *#
#*  You should have received a copy of the Apache-2.0 license                *#
#*  along with SoPlex; see the file LICENSE. If not email soplex@zib.de.     *#
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
#    HARDMEMLIMIT - hard memory limit given to slurm or the shell
#    HARDTIMELIMIT - hard time limit given to slurm
#    DEBUGTOOLCMD - debug tool command to prefix the executable with

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

# we add 10% to the hard memory limit and additional 100MB to the hard memory limit
HARDMEMLIMIT=$(((MEMLIMIT + 100) + (MEMLIMIT / 10)))

# in case of qsub queue the memory is measured in bytes
if test "${QUEUETYPE}" = "qsub"
then
    HARDMEMLIMIT=$((HARDMEMLIMIT * 1024000))
fi

# we add 100% to the hard time limit and additional 600 seconds in case of small time limits
HARDTIMELIMIT=$(((TIMELIMIT + 600) + TIMELIMIT))

# in case of srun the time needs to be formatted
if test "${QUEUETYPE}" != "qsub"
then
    MYMINUTES=0
    MYHOURS=0
    MYDAYS=0

    #calculate seconds, minutes, hours and days
    MYSECONDS=$((HARDTIMELIMIT % 60))
    TMP=$((HARDTIMELIMIT / 60))
    if test "${TMP}" != "0"
    then
        MYMINUTES=$((TMP % 60))
        TMP=$((TMP / 60))
        if test "${TMP}" != "0"
        then
            MYHOURS=$((TMP % 24))
            MYDAYS=$((TMP / 24))
        fi
    fi
    #format seconds to have two characters
    if test "${MYSECONDS}" -lt 10
    then
        MYSECONDS=0"${MYSECONDS}"
    fi
    #format minutes to have two characters
    if test "${MYMINUTES}" -lt 10
    then
        MYMINUTES=0"${MYMINUTES}"
    fi
    #format hours to have two characters
    if test "${MYHOURS}" -lt 10
    then
        MYHOURS=0"${MYHOURS}"
    fi
    #format HARDTIMELIMT
    if test "${MYDAYS}" = "0"
    then
        HARDTIMELIMIT="${MYHOURS}":"${MYMINUTES}":"${MYSECONDS}"
    else
        HARDTIMELIMIT="${MYDAYS}"-"${MYHOURS}":"${MYMINUTES}":"${MYSECONDS}"
    fi
fi

# check if the test run should be processed in a debug tool environment
if test "${DEBUGTOOL}" = "valgrind"
then
    DEBUGTOOLCMD="valgrind --log-fd=1 --leak-check=full --num-callers=20 --suppressions=${SOPLEXPATH}/../suppressions.valgrind "
elif test "${DEBUGTOOL}" = "rr"
then
    DEBUGTOOLCMD="rr record --chaos -o RRTRACEFOLDER_PLACEHOLDER.rr "
elif test "${DEBUGTOOL}" = "gdb"
then
    # leave a place holder for the error file to log to, which gets replaced in runcluster.sh
    DEBUGTOOLCMD='gdb -batch-silent -return-child-result -ex "run" -ex "set logging file ERRFILE_PLACEHOLDER" -ex "set logging on" -ex "thread apply all bt full" --args '
elif test "${DEBUGTOOL}" = "perf"
then
    DEBUGTOOLCMD="perf stat -e PERFEVENTS_PLACEHOLDER -x , -o PERFFILE_PLACEHOLDER -- "
else
    DEBUGTOOLCMD=""
fi

