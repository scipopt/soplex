#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#*                                                                           *#
#*                  This file is part of the class library                   *#
#*       SoPlex --- the Sequential object-oriented simPlex.                  *#
#*                                                                           *#
#*    Copyright (C) 1996-2014 Konrad-Zuse-Zentrum                            *#
#*                            fuer Informationstechnik Berlin                *#
#*                                                                           *#
#*  SoPlex is distributed under the terms of the ZIB Academic Licence.       *#
#*                                                                           *#
#*  You should have received a copy of the ZIB Academic License              *#
#*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  *#
#*                                                                           *#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#
# Call with "make testcluster"
#
# The queue is passed via $QUEUE (possibly defined in a local makefile in scip/make/local).
#
# For each run, we can specify the number of nodes reserved for a run via $PPN. If tests runs
# with valid time measurements should be executed, this number should be chosen in such a way
# that a job is run on a single computer, i.e., in general, $PPN should equal the number of cores
# of each computer. Of course, the value depends on the specific computer/queue.
#
# To get the result files call "./evalcheck_cluster.sh
# results/check.$TSTNAME.$BINNAME.$SETNAME.eval in directory check/
# This leads to result files
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.out
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.res
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.err

TSTNAME=$1
BINNAME=$2
BINID=$3
SETTINGS=$4
TIMELIMIT=$5
MEMLIMIT=$6
CONTINUE=$7
QUEUETYPE=$8
QUEUE=$9
PPN=${10}
CLIENTTMPDIR=${11}
NOWAITCLUSTER=${12}
EXCLUSIVE=${13}


# check whether all variables are defined
if test -z $EXCLUSIVE
then
    echo Skipping test since variable EXCLUSIVE is not defined.
    exit 1;
fi

# get current SOPLEX path
SOPLEXPATH=`pwd`

if test ! -e $SOPLEXPATH/results
then
    mkdir $SOPLEXPATH/results
fi

# check if binary exists
if test ! -e $SOPLEXPATH/../$BINNAME
then
    echo Skipping test since the binary $BINNAME does not exist.
    exit
fi

# check if queue has been defined
if test "$QUEUE" = ""
then
    echo Skipping test since the queue name has not been defined.
    exit
fi

# check if number of nodes has been defined
if test "$PPN" = ""
then
    echo Skipping test since the number of nodes has not been defined.
    exit
fi

# check if the slurm blades should be used exclusively
if test "$EXCLUSIVE" = "true"
then
    EXCLUSIVE=" --exclusive"
else
    EXCLUSIVE=""
fi

# we use the time limit as hard time limit because QSopt_ex's time limit does not work
HARDTIMELIMIT=$TIMELIMIT
HARDMEMLIMIT=$MEMLIMIT

# check whether there is enough memory on the host system, otherwise we need to submit from the target system
if test "$QUEUETYPE" = "srun"
then
    HOSTMEM=`ulimit -m`
    if test "$HOSTMEM" != "unlimited"
    then
        if [ `expr $HARDMEMLIMIT \* 1024` -gt $HOSTMEM ]
        then
            echo "Not enough memory on host system - please submit from target system (e.g. ssh opt201)."
            exit
        fi
    fi
fi

# in case of qsub queue the memory is measured in kB and in case of srun the time needs to be formatted
if test  "$QUEUETYPE" = "qsub"
then
    HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`
else
    MYMINUTES=0
    MYHOURS=0
    MYDAYS=0

    #calculate seconds, minutes, hours and days
    MYSECONDS=`expr $HARDTIMELIMIT % 60`
    TMP=`expr $HARDTIMELIMIT / 60`
    if test "$TMP" != "0"
    then
        MYMINUTES=`expr $TMP % 60`
        TMP=`expr $TMP / 60`
        if test "$TMP" != "0"
        then
            MYHOURS=`expr $TMP % 24`
            MYDAYS=`expr $TMP / 24`
        fi
    fi
    #format seconds to have two characters
    if test ${MYSECONDS} -lt 10
    then
        MYSECONDS=0${MYSECONDS}
    fi
    #format minutes to have two characters
    if test ${MYMINUTES} -lt 10
    then
        MYMINUTES=0${MYMINUTES}
    fi
    #format hours to have two characters
    if test ${MYHOURS} -lt 10
    then
        MYHOURS=0${MYHOURS}
    fi
    #format HARDTIMELIMT
    if test ${MYDAYS} = "0"
    then
        HARDTIMELIMIT=${MYHOURS}:${MYMINUTES}:${MYSECONDS}
    else
        HARDTIMELIMIT=${MYDAYS}-${MYHOURS}:${MYMINUTES}:${MYSECONDS}
    fi
fi

#define clusterqueue, which might not be the QUEUE, cause this might be an alias for a bunch of QUEUEs
CLUSTERQUEUE=$QUEUE

NICE=""
ACCOUNT="mip"

if test $CLUSTERQUEUE = "dbg"
then
    CLUSTERQUEUE="mip-dbg,telecom-dbg"
    ACCOUNT="mip-dbg"
elif test $CLUSTERQUEUE = "telecom-dbg"
then
    ACCOUNT="mip-dbg"
elif test $CLUSTERQUEUE = "mip-dbg"
then
    ACCOUNT="mip-dbg"
elif test $CLUSTERQUEUE = "opt-low"
then
    CLUSTERQUEUE="opt"
    NICE="--nice=10000"
fi


# counter to define file names for a test set uniquely
COUNT=0

for i in `cat $SOPLEXPATH/testset/$TSTNAME.test`
do
    if test "$i" = "DONE"
        then
        break
    fi

    # increase the index for the inctance tried to solve, even if the filename does not exist
    COUNT=`expr $COUNT + 1`

    # check if problem instance exists
    if test -f $SOPLEXPATH/$i
    then

        # the cluster queue has an upper bound of 2000 jobs; if this limit is
        # reached the submitted jobs are dumped; to avoid that we check the total
        # load of the cluster and wait until it is save (total load not more than
        # 1900 jobs) to submit the next job.
        if test "$NOWAITCLUSTER" != "1"
        then
            if test  "$QUEUETYPE" != "qsub"
            then
                echo "waitcluster does not work on slurm cluster"
            fi
            ./waitcluster.sh 1600 $QUEUE 200
        fi

        SHORTFILENAME=`basename $i .gz`
        SHORTFILENAME=`basename $SHORTFILENAME .mps`
        SHORTFILENAME=`basename $SHORTFILENAME .lp`
        SHORTFILENAME=`basename $SHORTFILENAME .opb`

        FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTFILENAME.$BINID.$QUEUE.$SETTINGS
        BASENAME=$SOPLEXPATH/results/$FILENAME

        TMPFILE=$BASENAME.tmp
        SETFILE=$BASENAME.set

        # in case we want to continue we check if the job was already performed
        if test "$CONTINUE" != "false"
        then
            if test -e results/$FILENAME.out
            then
                echo skipping file $i due to existing output file $FILENAME.out
                continue
            fi
        fi

        # additional environment variables needed by runcluster.sh
        export SOLVERPATH=$SOPLEXPATH
        export EXECNAME=$SOPLEXPATH/../$BINNAME
        export BASENAME=$FILENAME
        export FILENAME=$i
        export TIMELIMIT=$TIMELIMIT
        export SETTINGS=$SETTINGS
        export INSTANCE=$SOPLEXPATH/$i
        export CLIENTTMPDIR=$CLIENTTMPDIR

        # check queue type
        if test  "$QUEUETYPE" = "srun"
        then
            sbatch --job-name=QSO-$SHORTFILENAME --mem=$HARDMEMLIMIT -p $CLUSTERQUEUE -A $ACCOUNT --time=${HARDTIMELIMIT} ${NICE} ${EXCLUSIVE} --output=/dev/null runcluster_qsoptex.sh
        else
            # -V to copy all environment variables
            qsub -l walltime=$HARDTIMELIMIT -l mem=$HARDMEMLIMIT -l nodes=1:ppn=$PPN -N SOPLEX$SHORTFILENAME -V -q $CLUSTERQUEUE -o /dev/null -e /dev/null runcluster_qsoptex.sh
        fi
    else
        echo "input file "$SOPLEXPATH/$i" not found!"
    fi
done
