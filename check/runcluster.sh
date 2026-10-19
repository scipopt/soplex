#!/usr/bin/env bash
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

# Execute binary and write logfiles, gets called by check_cluster.sh or run_group.sh.
#
# When called with a phase argument (setup/solve), only that phase runs.
# Without argument (or "all"), all phases run in sequence.

PHASE="${1:-all}"

# check if tmp-path exists
if test ! -d "${CLIENTTMPDIR}"
then
    echo Skipping test since the path for the tmp-dir does not exist.
    exit
fi

OUTFILE="${CLIENTTMPDIR}/${BASENAME}.out"
ERRFILE="${CLIENTTMPDIR}/${BASENAME}.err"
BASFILE="${CLIENTTMPDIR}/${BASENAME}.bas"
PERFFILE="${CLIENTTMPDIR}/${BASENAME}.perf"
SETFILE="${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.set"

# function to move the result files to ${OUTPUTDIR}
function cleanup {
    if [[ "${PHASE}" == "all" || "${PHASE}" == "solve" ]]
    then
        mv "${OUTFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.out"
        mv "${ERRFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.err"
        if test -f "${BASFILE}"
        then
            mv "${BASFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.bas"
        fi
        if test -f "${PERFFILE}"
        then
            mv "${PERFFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.perf"
        fi
    fi
}

# ensure the result files are moved to ${OUTPUTDIR} when exiting (normally or due to abort/interrupt)
trap cleanup EXIT

if [[ "${PHASE}" == "all" || "${PHASE}" == "setup" ]]
then
    uname -a                               > "${OUTFILE}"
    uname -a                               > "${ERRFILE}"

    echo "@01 ${FILENAME}"                >> "${OUTFILE}"
    echo "@01 ${FILENAME}"                >> "${ERRFILE}"
    echo "-----------------------------"  >> "${OUTFILE}"
    date                                  >> "${OUTFILE}"
    date                                  >> "${ERRFILE}"
    echo "-----------------------------"  >> "${OUTFILE}"
    date +"@03 %s"                        >> "${OUTFILE}"
    echo "@05 ${TIMELIMIT}"               >> "${OUTFILE}"

    # system state before solve: available memory (kB), cumulative CPU ticks (idle total)
    if test -r /proc/meminfo && test -r /proc/stat
    then
        echo "@07 $(awk '/MemAvailable/ {print $2}' /proc/meminfo) $(awk '/^cpu / {printf "%d %d\n", $5, $2+$3+$4+$5+$6+$7+$8+$9}' /proc/stat)" >> "${OUTFILE}"
    fi

fi

# resolve placeholders in the executable name (done by the caller in setup/solve mode)
if [[ "${PHASE}" == "all" ]]
then
    . ./prepare_execname.sh
fi

if [[ "${PHASE}" == "all" || "${PHASE}" == "solve" ]]
then
    SOLVECMD="${EXECNAME} --saveset=\"${SETFILE}\" --writebas=\"${BASFILE}\" \"${INSTANCE}\""

    # hold the core until all steps of the group are submitted
    if test -n "${HOLD}"
    then
        SOLVECMD="bash -c 'sleep ${HOLD}; ${SOLVECMD}'"
    fi

    eval "${SRUN:+${SRUN} }${SOLVECMD}" >> "${OUTFILE}" 2>> "${ERRFILE}"

    # system state after solve: same format as @07
    if test -r /proc/meminfo && test -r /proc/stat
    then
        echo "@08 $(awk '/MemAvailable/ {print $2}' /proc/meminfo) $(awk '/^cpu / {printf "%d %d\n", $5, $2+$3+$4+$5+$6+$7+$8+$9}' /proc/stat)" >> "${OUTFILE}"
    fi
    echo "@09 ${TARGETFREQ:-0}"           >> "${OUTFILE}"
    date +"@04 %s"                        >> "${OUTFILE}"
    echo "-----------------------------"  >> "${OUTFILE}"
    date                                  >> "${OUTFILE}"
    echo "-----------------------------"  >> "${OUTFILE}"
    date                                  >> "${ERRFILE}"
fi
