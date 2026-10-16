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

# launches multiple solver instances in parallel on one exclusive node.
# used by check_cluster.sh when EXCLUSIVE=auto to pack instances onto one node.
# setup runs first unbound, then solve steps are staggered for stable core binding.
#
# usage: run_group.sh SRUN_PREFIX NUM_INSTANCES [FIELDS...]
#
# FIELDS contains NUM_INSTANCES groups of 10 positional arguments each:
#   SOLVERPATH EXECNAME BASENAME FILENAME OUTPUTDIR TIMELIMIT SETTINGS INSTANCE CLIENTTMPDIR TARGETFREQ

SRUN_PREFIX="$1"
NUM_INSTANCES="$2"
shift 2

# encode per-instance parameters with prepared EXECNAME
for ((i = 0; i < NUM_INSTANCES; i++))
do
    EXECNAME="$2"
    ERRFILE="$9/$3.err"
    PERFFILE="$9/$3.perf"
    . ./prepare_execname.sh

    printf -v "INST_${i}" 'export SOLVERPATH=%q EXECNAME=%q BASENAME=%q FILENAME=%q OUTPUTDIR=%q TIMELIMIT=%q SETTINGS=%q INSTANCE=%q CLIENTTMPDIR=%q TARGETFREQ=%q' \
        "$1" "${EXECNAME}" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}"
    export "INST_${i}"
    shift 10
done

# phase 1: setup solver environments
for ((i = 0; i < NUM_INSTANCES; i++))
do
    (n=INST_${i}; eval "${!n}"; ./runcluster.sh setup) &
done
wait

# hold cores until all steps submitted
TICK=0.03125
HOLD=$(awk "BEGIN {printf \"%.5f\", ${NUM_INSTANCES} * ${TICK}}")

# phase 2: solve and cleanup
for ((i = 0; i < NUM_INSTANCES; i++))
do
    sleep ${TICK}
    (n=INST_${i}; eval "${!n}"; export SRUN="${SRUN_PREFIX}"; export HOLD="${HOLD}"; ./runcluster.sh solve) &
done
wait
