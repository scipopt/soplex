#!/bin/bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the class library                   *
#*       SoPlex --- the Sequential object-oriented simPlex.                  *
#*                                                                           *
#*  Copyright (c) 1996-2025 Zuse Institute Berlin (ZIB)                      *
#*                                                                           *
#*  Licensed under the Apache License, Version 2.0 (the "License");          *
#*  you may not use this file except in compliance with the License.         *
#*  You may obtain a copy of the License at                                  *
#*                                                                           *
#*      http://www.apache.org/licenses/LICENSE-2.0                           *
#*                                                                           *
#*  Unless required by applicable law or agreed to in writing, software      *
#*  distributed under the License is distributed on an "AS IS" BASIS,        *
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
#*  See the License for the specific language governing permissions and      *
#*  limitations under the License.                                           *
#*                                                                           *
#*  You should have received a copy of the Apache-2.0 license                *
#*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de.  *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# configures environment variables for cluster runs.
# to be invoked inside a check_cluster.sh script
# This script cancels the process if required variables are not correctly set

# input variables - should be passed to this script
# QUEUE=${1}    # the name of the cluster (M620, dbg, telecom-dbg, mip-dbg, opt-low, opt)
# PPN=${2}      # number of cluster nodes to use
# EXCLUSIVE=${3} # should cluster nodes be blocked for other users while the jobs are running?
# QUEUETYPE=${4} # either 'srun' or 'qsub'

# new environment variables defined by this script:
#   HARDTIMELIMIT - hard time limit
#   NICE          - cluster job nice flag
#   ACCOUNT       - cluster job account flag
#   CLUSTERQUEUE  - list of cluster queues

NICE=""
if [[ "$(uname -n)" =~ htc ]]; then
  # z1 cluster
  ACCOUNT="optimi_integer"
else
  # opt machines
  ACCOUNT="mip"
fi
CLUSTERQUEUE="${QUEUE}"

# check if queue has been defined
if test "${QUEUE}" = ""
then
    echo "Skipping test since the queue name has not been defined."
    exit 1
fi

# check if number of nodes has been defined
if test "${PPN}" = ""
then
    echo "Skipping test since the number of nodes has not been defined."
    exit 1
fi

# we add 100% to the hard time limit and additional 600 seconds in case of small time limits
# NOTE: the jobs should have a hard running time of more than 5 minutes; if not so, these
#       jobs get automatically assigned in the "exrpess" queue; this queue has only 4 CPUs
#       available
HARDTIMELIMIT=$(((TIMELIMIT + 600) + TIMELIMIT))

# we add 10% to the hard memory limit and additional 100MB to the hard memory limit
HARDMEMLIMIT=$(((MEMLIMIT + 100) + (MEMLIMIT / 10)))


# check whether there is enough memory on the host system, otherwise we need to submit from the target system
if test "${QUEUETYPE}" = "srun"
then
    HOSTMEM=$(ulimit -m)
    if test "${HOSTMEM}" != "unlimited"
    then
        if [ $((HARDMEMLIMIT * 1024)) -gt "${HOSTMEM}" ]
        then
            echo "Not enough memory on host system - please submit from target system (e.g. ssh opt201)."
            exit 1
        fi
    fi
fi

# in case of qsub queue the memory is measured in kB and in case of srun the time needs to be formatted
if test  "${QUEUETYPE}" = "qsub"
then
    HARDMEMLIMIT=$((HARDMEMLIMIT * 1024000))
else
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


#define clusterqueue, which might not be the QUEUE, because this might be an alias for a bunch of QUEUEs
if test "${CLUSTERQUEUE}" = "dbg"
then
    CLUSTERQUEUE="mip-dbg"
    ACCOUNT="mip-dbg"
elif test "${CLUSTERQUEUE}" = "mip-dbg"
then
    ACCOUNT="mip-dbg"
elif test "${CLUSTERQUEUE}" = "opt-low"
then
    CLUSTERQUEUE="opt"
    NICE="--nice=10000"

    # wakeup the cluster
    make --makefile=wakeup-slurm wake_opt
elif test "${CLUSTERQUEUE}" = "M620-low"
then
    NICE="--nice=10000"
    CLUSTERQUEUE="M620"

    # wakeup the cluster
    make --makefile=wakeup-slurm wake_M620
elif test "${CLUSTERQUEUE}" = "M620v3-low"
then
    NICE="--nice=10000"
    CLUSTERQUEUE="M620v3"

    # wakeup the cluster
    make --makefile=wakeup-slurm wake_M620v3
elif test "${CLUSTERQUEUE}" = "M630-low"
then
    NICE="--nice=10000"
    CLUSTERQUEUE="M630"

    # wakeup the cluster
    make --makefile=wakeup-slurm wake_M630
elif test "${CLUSTERQUEUE}" = "M620x"
then
    CLUSTERQUEUE="M620,M620v2,M620v3"

    # wakeup the cluster
    make --makefile=wakeup-slurm wake_M620
    make --makefile=wakeup-slurm wake_M620v2
    make --makefile=wakeup-slurm wake_M620v3
elif test "${CLUSTERQUEUE}" = "M640-low"
then
    NICE="--nice=10000"
    CLUSTERQUEUE="M640"

    # wakeup the cluster
    make --makefile=wakeup-slurm wake_M640
elif test "${CLUSTERQUEUE}" = "moskito"
then
    ACCOUNT="dopt"
fi

# check if the slurm blades should be used exclusively
if test "${EXCLUSIVE}" = "true"
then
    EXCLUSIVE=" --exclusive"
else
    EXCLUSIVE=""
fi


CONSTRAINT=""
if test "${CLUSTERQUEUE}" = "Gold6338"
then
    CONSTRAINT="Gold6338"
    CLUSTERQUEUE="big"
elif test "${CLUSTERQUEUE}" = "Gold6342"
then
    CONSTRAINT="Gold6342"
    CLUSTERQUEUE="big"
elif test "${CLUSTERQUEUE}" = "M640v2"
then
    CONSTRAINT="Gold5222"
    CLUSTERQUEUE="opt_int"
elif test "${CLUSTERQUEUE}" = "M640"
then
    CONSTRAINT="Gold5122"
    CLUSTERQUEUE="opt_int"
fi
