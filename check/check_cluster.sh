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
# Call with "make testcluster"
#
# The queue is passed via ${QUEUE} (possibly defined in a local makefile in soplex/make/local).
#
# For each run, we can specify the number of nodes reserved for a run via ${PPN}. With EXCLUSIVE=auto
# the number of instances packed onto one exclusive node is computed automatically, so that each
# instance is bound to its own physical core for valid time measurements.
#
# To get the result files call in directory check/
# ./evalcheck_cluster.sh ${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${SETNAME}.eval
# This leads to result files
#  - ${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${SETNAME}.out
#  - ${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${SETNAME}.res
#  - ${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${SETNAME}.err
#
# To get verbose output from Slurm, have SRUN_FLAGS="-v -v" set in your environment.

TSTNAME="${1}"
EXECUTABLE="${2}"
SETTINGS="${3}"
BINID="${4}"
OUTPUTDIR="${5}"
TIMELIMIT="${6}"
MEMLIMIT="${7}"
CONTINUE="${8}"
QUEUETYPE="${9}"
QUEUE="${10}"
PPN="${11}"
CLIENTTMPDIR="${12}"
NOWAITCLUSTER="${13}"
EXCLUSIVE="${14}"
SEEDS="${15}"
DEBUGTOOL="${16}"
CLUSTERNODES="${17}"
EXCLUDENODES="${18}"
SLURMACCOUNT="${19}"

# check if all variables defined (by checking the last one)
if test -z "${SLURMACCOUNT}"
then
    echo Skipping test since not all variables are defined
    echo "TSTNAME       = ${TSTNAME}"
    echo "EXECUTABLE    = ${EXECUTABLE}"
    echo "SETTINGS      = ${SETTINGS}"
    echo "BINID         = ${BINID}"
    echo "OUTPUTDIR     = ${OUTPUTDIR}"
    echo "TIMELIMIT     = ${TIMELIMIT}"
    echo "MEMLIMIT      = ${MEMLIMIT}"
    echo "CONTINUE      = ${CONTINUE}"
    echo "QUEUETYPE     = ${QUEUETYPE}"
    echo "QUEUE         = ${QUEUE}"
    echo "PPN           = ${PPN}"
    echo "CLIENTTMPDIR  = ${CLIENTTMPDIR}"
    echo "NOWAITCLUSTER = ${NOWAITCLUSTER}"
    echo "EXCLUSIVE     = ${EXCLUSIVE}"
    echo "SEEDS         = ${SEEDS}"
    echo "DEBUGTOOL     = ${DEBUGTOOL}"
    echo "CLUSTERNODES  = ${CLUSTERNODES}"
    echo "EXCLUDENODES  = ${EXCLUDENODES}"
    echo "SLURMACCOUNT  = ${SLURMACCOUNT}"
    exit 1
fi

# call routines for creating the result directory, checking for existence
# of passed settings, etc
# defines the following environment variables: SOPLEXPATH, EXECUTABLE, FULLTSTNAME, SOLUFILE, SETTINGSFILE, HARDMEMLIMIT, HARDTIMELIMIT, DEBUGTOOLCMD
. ./configuration_set.sh

# create output directory
mkdir -p "${SOPLEXPATH}/${OUTPUTDIR}"

# check if binary exists
if test ! -e "${EXECUTABLE}"
then
    echo "Skipping test since the binary ${EXECUTABLE} does not exist."
    exit 1
fi

# configure cluster-related environment variables
# defines the following environment variables: NICE, ACCOUNT, CLUSTERQUEUE, CONSTRAINT, CPUFREQ, TARGETFREQ, EXCLUSIVE, AUTO_PPN_PENDING
. ./configuration_cluster.sh

# prefix the executable with the debug tool command (placeholders resolved in runcluster.sh)
EXECNAME="${DEBUGTOOLCMD}${EXECUTABLE}"

# resolve the default slurm account to the one chosen by configuration_cluster.sh
if test "${SLURMACCOUNT}" = "default"
then
    SLURMACCOUNT="${ACCOUNT}"
fi

# compute node filter flags (CLUSTERNODES/EXCLUDENODES don't change per instance)
if test "${CLUSTERNODES}" = "all" && test "${EXCLUDENODES}" = "none"
then
    NODE_FLAGS=""
elif test "${CLUSTERNODES}" != "all" && test "${EXCLUDENODES}" = "none"
then
    NODE_FLAGS="-w ${CLUSTERNODES}"
elif test "${CLUSTERNODES}" = "all" && test "${EXCLUDENODES}" != "none"
then
    NODE_FLAGS="-x ${EXCLUDENODES}"
else
    NODE_FLAGS="-w ${CLUSTERNODES} -x ${EXCLUDENODES}"
fi

# query node topology to compute instances-per-node
if test "${AUTO_PPN_PENDING}" -eq 1
then
    SINFO_CMD="sinfo --noheader -N -p ${CLUSTERQUEUE} -O"
    if test -n "${CONSTRAINT}"
    then
        SINFO_FILTER="grep ${CONSTRAINT}"
    else
        SINFO_FILTER="cat"
    fi
    NODE_SOCKETS=$(${SINFO_CMD} sockets,features | ${SINFO_FILTER} | awk '{print $1}' | sort -n | head -1 | tr -d ' ')
    NODE_CORES=$(${SINFO_CMD} cores,features | ${SINFO_FILTER} | awk '{print $1}' | sort -n | head -1 | tr -d ' ')

    if test "${NODE_SOCKETS:-0}" -gt 0 && test "${NODE_CORES:-0}" -gt 0
    then
        NODE_CPUS=$(( NODE_CORES * NODE_SOCKETS ))
        PPN_CPU="${NODE_CPUS}"
    else
        echo "Warning: could not determine node cores for partition ${CLUSTERQUEUE}."
    fi
    NODE_MEM=$(${SINFO_CMD} memory,features | ${SINFO_FILTER} | awk '{print $1}' | sort -n | head -1 | tr -d ' ')
    if test "${NODE_MEM:-0}" -gt 0
    then
        if test "${NODE_MEM}" -lt "${HARDMEMLIMIT}"
        then
            echo "Error: memory exceeds available node storage resource: NODE_MEM=${NODE_MEM}, HARDMEMLIMIT=${HARDMEMLIMIT}."
            exit 1
        fi
        PPN_MEM=$(( NODE_MEM / HARDMEMLIMIT ))
    else
        echo "Warning: could not determine node memory for partition ${CLUSTERQUEUE}."
    fi
    if test "${PPN}" -gt 0
    then
        if test -n "${PPN_CPU}" && test "${PPN}" -gt "${PPN_CPU}"
        then
            echo "Error: PPN=${PPN} exceeds PPN_CPU=${PPN_CPU} (NODE_CPUS=${NODE_CPUS})."
            exit 1
        fi
        if test -n "${PPN_MEM}" && test "${PPN}" -gt "${PPN_MEM}"
        then
            echo "Error: PPN=${PPN} exceeds PPN_MEM=${PPN_MEM} (NODE_MEM=${NODE_MEM}, HARDMEMLIMIT=${HARDMEMLIMIT})."
            exit 1
        fi
        AUTO_PPN="${PPN}"
    else
        AUTO_PPN="${PPN_CPU}"
        if test -z "${AUTO_PPN}" || test -n "${PPN_MEM}" && test "${PPN_MEM}" -lt "${AUTO_PPN}"
        then
            AUTO_PPN="${PPN_MEM}"
        fi
        if test -z "${AUTO_PPN}"
        then
            AUTO_PPN=1
        fi
    fi

    # srun prefix for the within-job steps
    if test "${CLUSTERQUEUE}" != "moskito" && test "${CLUSTERQUEUE}" != "prio"
    then
        AUTO_SRUN="srun --exact -n 1 -c 1 --mem=${HARDMEMLIMIT} --cpu_bind=verbose,cores --propagate=STACK ${SRUN_FLAGS}"
    fi
fi

# non-auto srun uses one instance per sbatch
if test "${QUEUETYPE}" = "srun" && test "${AUTO_PPN_PENDING}" -eq 0 && test "${PPN}" -gt 1
then
    echo "Error: PPN=${PPN} requires EXCLUSIVE=auto."
    exit 1
fi

# auto mode batch accumulator
AUTO_BATCH_COUNT=0          # number of instances accumulated so far
AUTO_BATCH_JOBNAME=""       # job name of first instance in group
AUTO_BATCH_NAMES=""         # space-separated instance names for echo
AUTO_BATCH_ARGS=()          # positional arguments for run_group.sh (10 per instance)

# flush the current auto mode batch as a single exclusive sbatch job
flush_auto_batch() {
    if test "${AUTO_BATCH_COUNT}" -eq 0
    then
        return
    fi
    GROUP_MEM=$(( AUTO_BATCH_COUNT * HARDMEMLIMIT ))
    echo sbatch --job-name="${AUTO_BATCH_JOBNAME}" --constraint="${CONSTRAINT}" --mem="${GROUP_MEM}" -n "${AUTO_BATCH_COUNT}" -c 1 --distribution=cyclic:cyclic -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${HARDTIMELIMIT}" --cpu-freq=${CPUFREQ} --exclusive ${NODE_FLAGS} --output=/dev/null run_group.sh
    echo "instances:${AUTO_BATCH_NAMES}"
    sbatch --job-name="${AUTO_BATCH_JOBNAME}" --constraint="${CONSTRAINT}" --mem="${GROUP_MEM}" -n "${AUTO_BATCH_COUNT}" -c 1 --distribution=cyclic:cyclic -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${HARDTIMELIMIT}" --cpu-freq=${CPUFREQ} --exclusive ${NODE_FLAGS} --output=/dev/null run_group.sh "${AUTO_SRUN}" "${AUTO_BATCH_COUNT}" "${AUTO_BATCH_ARGS[@]}"
    # reset accumulator
    AUTO_BATCH_COUNT=0
    AUTO_BATCH_JOBNAME=""
    AUTO_BATCH_NAMES=""
    AUTO_BATCH_ARGS=()
}

# run different random seeds
for (( s=0; s<=${SEEDS}; s++ ))
do
    EVALFILE="${SOPLEXPATH}/${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${QUEUE}.${SETTINGS}"
    if (( s > 0 )); then
        EVALFILE="${EVALFILE}-s${s}"
    fi
    EVALFILE="${EVALFILE}.eval"
    echo > "${EVALFILE}"

    METAFILE="${EVALFILE%.eval}.meta"
    echo "@Seed ${s}"              >  "${METAFILE}"
    echo "@Settings ${SETTINGS}"   >> "${METAFILE}"
    echo "@TstName ${TSTNAME}"     >> "${METAFILE}"
    echo "@BinName ${BINID}"       >> "${METAFILE}"
    echo "@MemLimit ${MEMLIMIT}"   >> "${METAFILE}"
    echo "@TimeLimit ${TIMELIMIT}" >> "${METAFILE}"
    echo "@Queue ${QUEUE}"         >> "${METAFILE}"
    echo "@Exclusive ${EXCLUSIVE}" >> "${METAFILE}"

    # counter to define file names for a test set uniquely
    COUNT=0

    # iterate over instances in .test file
    for i in $(cat "${FULLTSTNAME}")
    do
        if test "${i}" = "DONE"
            then
            break
        fi

        # increase the index for the instance tried to solve, even if the filename does not exist
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

            SHORTNAME="${USER}.${TSTNAME}.${COUNT}_${SHORTFILENAME}.${BINID}.${QUEUE}.${SETTINGS}"
            if (( s > 0 )); then
                SHORTNAME="${SHORTNAME}-s${s}"
            fi
            RESBASE="${SOPLEXPATH}/${OUTPUTDIR}/${SHORTNAME}"

            echo "${RESBASE}" >> "${EVALFILE}"

            # in case we want to continue we check if the job was already performed
            if test "${CONTINUE}" != "false"
            then
                if test -e "${OUTPUTDIR}/${SHORTNAME}.out"
                then
                    echo "skipping file ${i} due to existing output file ${SHORTNAME}.out"
                    continue
                fi
            fi

            # additional environment variables needed by runcluster.sh
            export SOLVERPATH="${SOPLEXPATH}"
            export EXECNAME
            export BASENAME="${SHORTNAME}" # careful, runcluster.sh uses BASENAME in a different sense
            export FILENAME="${i}"         # careful, runcluster.sh uses FILENAME in a different sense
            export OUTPUTDIR="${OUTPUTDIR}"
            export TIMELIMIT="${TIMELIMIT}"
            export SETTINGS="${SETTINGS}"
            export INSTANCE="${SOPLEXPATH}/${i}"
            export CLIENTTMPDIR="${CLIENTTMPDIR}"
            export TARGETFREQ

            # Create testset, without printing output
            "${EXECUTABLE}" --loadset="${SETTINGSFILE}" -t"${TIMELIMIT}" --saveset="${SOLVERPATH}/${OUTPUTDIR}/${SHORTNAME}.set" > /dev/null

            if test  "${QUEUETYPE}" = "srun"
            then
                if test "${AUTO_PPN_PENDING}" -eq 1
                then
                    # auto mode: accumulate instance arguments into current batch
                    AUTO_BATCH_ARGS+=("${SOLVERPATH}" "${EXECNAME}" "${SHORTNAME}" "${i}" "${OUTPUTDIR}" "${TIMELIMIT}" "${SETTINGS}" "${SOLVERPATH}/${i}" "${CLIENTTMPDIR}" "${TARGETFREQ}")
                    AUTO_BATCH_COUNT=$(( AUTO_BATCH_COUNT + 1 ))
                    AUTO_BATCH_NAMES="${AUTO_BATCH_NAMES} ${SHORTFILENAME}"
                    if test -z "${AUTO_BATCH_JOBNAME}"
                    then
                        AUTO_BATCH_JOBNAME="SPX-${SHORTFILENAME}"
                    fi
                    # flush when batch is full
                    if test "${AUTO_BATCH_COUNT}" -ge "${AUTO_PPN}"
                    then
                        flush_auto_batch
                    fi
                else
                    if test "${CLUSTERQUEUE}" != "moskito" && test "${CLUSTERQUEUE}" != "prio"
                    then
                        export SRUN="srun --exact -n 1 -c 1 --mem=${HARDMEMLIMIT} --propagate=STACK --cpu_bind=verbose,cores ${SRUN_FLAGS}"
                    fi
                    echo sbatch --job-name=SPX-"${SHORTFILENAME}" --constraint="${CONSTRAINT}" -n 1 -c 1 --mem="${HARDMEMLIMIT}" -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${HARDTIMELIMIT}" --cpu-freq=${CPUFREQ} ${EXCLUSIVE} ${NODE_FLAGS} --output=/dev/null runcluster.sh
                    sbatch --job-name=SPX-"${SHORTFILENAME}" --constraint="${CONSTRAINT}" -n 1 -c 1 --mem="${HARDMEMLIMIT}" -p "${CLUSTERQUEUE}" -A "${SLURMACCOUNT}" ${NICE} --time="${HARDTIMELIMIT}" --cpu-freq=${CPUFREQ} ${EXCLUSIVE} ${NODE_FLAGS} --output=/dev/null runcluster.sh
                fi
            else
                # -V to copy all environment variables
                qsub -l walltime="${HARDTIMELIMIT}" -l mem="${HARDMEMLIMIT}" -l nodes=1:ppn="${PPN}" -N SOPLEX"${SHORTFILENAME}" -V -q "${CLUSTERQUEUE}" -o /dev/null -e /dev/null runcluster.sh
            fi
        else
            echo "input file ${SOPLEXPATH}/${i} not found!"
        fi
    done

    # flush any remaining instances in auto mode at end of each seed
    if test "${AUTO_PPN_PENDING}" -eq 1 && test "${AUTO_BATCH_COUNT}" -gt 0
    then
        flush_auto_batch
    fi
done
