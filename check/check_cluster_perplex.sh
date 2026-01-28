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
#
# Call with "make testclusterperplex"
#
# The queue is passed via ${QUEUE} (possibly defined in a local makefile in scip/make/local).
#
# For each run, we can specify the number of nodes reserved for a run via ${PPN}. If tests runs
# with valid time measurements should be executed, this number should be chosen in such a way
# that a job is run on a single computer, i.e., in general, ${PPN} should equal the number of cores
# of each computer. Of course, the value depends on the specific computer/queue.
#
# To get the result files call in directory check/
# ./evalcheck_cluster.sh results/check.${TSTNAME}.${BINID}.${SETNAME}.eval
# This leads to result files
#  - results/check.${TSTNAME}.${BINID}.${SETNAME}.out
#  - results/check.${TSTNAME}.${BINID}.${SETNAME}.res
#  - results/check.${TSTNAME}.${BINID}.${SETNAME}.err

TSTNAME="${1}"
EXECUTABLE="${2}"
BINID="${3}"
SETTINGS="${4}"
TIMELIMIT="${5}"
MEMLIMIT="${6}"
CONTINUE="${7}"
QUEUETYPE="${8}"
QUEUE="${9}"
PPN="${10}"
CLIENTTMPDIR="${11}"
NOWAITCLUSTER="${12}"
EXCLUSIVE="${13}"
OUTPUTDIR="${14}"

# check if all variables defined (by checking the last one)
if test -z "${OUTPUTDIR}"
then
    echo Skipping test since not all variables are defined
    echo "TSTNAME       = ${TSTNAME}"
    echo "EXECUTABLE    = ${EXECUTABLE}"
    echo "BINID         = ${BINID}"
    echo "SETTINGS      = ${SETTINGS}"
    echo "TIMELIMIT     = ${TIMELIMIT}"
    echo "MEMLIMIT      = ${MEMLIMIT}"
    echo "CONTINUE      = ${CONTINUE}"
    echo "QUEUETYPE     = ${QUEUETYPE}"
    echo "QUEUE         = ${QUEUE}"
    echo "PPN           = ${PPN}"
    echo "CLIENTTMPDIR  = ${CLIENTTMPDIR}"
    echo "NOWAITCLUSTER = ${NOWAITCLUSTER}"
    echo "EXCLUSIVE     = ${EXCLUSIVE}"
    echo "OUTPUTDIR     = ${OUTPUTDIR}"
    exit 1;
fi

# call routines for creating the result directory, checking for existence
# of passed settings, etc
# defines the following environment variables: SOPLEXPATH, EXECUTABLE, FULLTSTNAME, SOLUFILE, SETTINGSFILE
. ./configuration_set.sh

# create output directory
mkdir -p "${SOPLEXPATH}/${OUTPUTDIR}"

# check if binary exists
if test ! -e "${EXECUTABLE}"
then
    echo "Skipping test since the binary ${EXECUTABLE} does not exist."
    exit
fi

# we use the time limit as hard time limit because perPlex does not accept a time limit
HARDTIMELIMIT="${TIMELIMIT}"
HARDMEMLIMIT="${MEMLIMIT}"

# configure cluster-related environment variables
# defines the following environment variables: NICE, ACCOUNT, CLUSTERQUEUE, HARDTIMELIMIT
. ./configuration_cluster.sh

# counter to define file names for a test set uniquely
COUNT=0

# iterate over instances in .test file
for i in $(cat "${FULLTSTNAME}")
do
    if test "${i}" = "DONE"
        then
        break
    fi

    # increase the index for the inctance tried to solve, even if the filename does not exist
    COUNT=$((COUNT + 1))

    # check if problem instance exists
    if test -f "${SOPLEXPATH}/${i}"
    then

        # the cluster queue has an upper bound of 2000 jobs; if this limit is
        # reached the submitted jobs are dumped; to avoid that we check the total
        # load of the cluster and wait until it is save (total load not more than
        # 1900 jobs) to submit the next job.
        if test "${NOWAITCLUSTER}" != "1"
        then
            if test  "${QUEUETYPE}" != "qsub"
            then
                echo "waitcluster does not work on slurm cluster"
            fi
            ./waitcluster.sh 1600 "${QUEUE}" 200
        fi

        SHORTFILENAME=$(basename "${i}" .gz)
        SHORTFILENAME=$(basename "${SHORTFILENAME}" .mps)
        SHORTFILENAME=$(basename "${SHORTFILENAME}" .lp)
        SHORTFILENAME=$(basename "${SHORTFILENAME}" .opb)

        FILENAME="${USER}.${TSTNAME}.${COUNT}_${SHORTFILENAME}.${BINID}.${QUEUE}.${SETTINGS}"
        BASENAME="${SOPLEXPATH}/results/${FILENAME}"

        TMPFILE="${BASENAME}.tmp"
        SETFILE="${BASENAME}.set"

        # in case we want to continue we check if the job was already performed
        if test "${CONTINUE}" != "false"
        then
            if test -e "results/${FILENAME}.out"
            then
                echo "skipping file ${i} due to existing output file ${FILENAME}.out"
                continue
            fi
        fi

        # additional environment variables needed by runcluster.sh
        export SOLVERPATH="${SOPLEXPATH}"
        export EXECNAME="${EXECUTABLE}"
        export BASENAME="${FILENAME}"
        export FILENAME="${i}"
        export TIMELIMIT="${TIMELIMIT}"
        export SETTINGS="${SETTINGS}"
        export INSTANCE="${SOPLEXPATH}/${i}"
        export CLIENTTMPDIR="${CLIENTTMPDIR}"

        # check queue type
        if test  "${QUEUETYPE}" = "srun"
        then
            sbatch --job-name=PPX-"${SHORTFILENAME}" --mem="${HARDMEMLIMIT}" -p "${CLUSTERQUEUE}" -A "${ACCOUNT}" --time="${HARDTIMELIMIT}" ${NICE} ${EXCLUSIVE} --output=/dev/null runcluster_perplex.sh
        else
            # -V to copy all environment variables
            qsub -l walltime="${HARDTIMELIMIT}" -l mem="${HARDMEMLIMIT}" -l nodes=1:ppn="${PPN}" -N SOPLEX"${SHORTFILENAME}" -V -q "${CLUSTERQUEUE}" -o /dev/null -e /dev/null runcluster_perplex.sh
        fi
    else
        echo "input file ${SOPLEXPATH}/${i} not found!"
    fi
done
