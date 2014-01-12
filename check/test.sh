#!/bin/sh

# solve a given testset with given settings and time limit
# parameters:
# 1: name of testset (has to be in check/testset)
# 2: path to soplex executable
# 3: name of settings (has to be in settings)
# 4: time limit

TEST=$1
TESTSET=testset/$TEST.test
SOLUNAME=testset/$TEST.solu

BINFILE=$2
BINNAME=`basename $2`

SETTINGS=$3
SETTINGSFILE=../settings/$SETTINGS.set

TIME=$4

OUTFILE=results/check.$TEST.$BINNAME.$SETTINGS.out
ERRFILE=results/check.$TEST.$BINNAME.$SETTINGS.err
RESFILE=results/check.$TEST.$BINNAME.$SETTINGS.res
SETFILE=results/check.$TEST.$BINNAME.$SETTINGS.set

# Create results directory
mkdir -p results

# Abort if files are missing
if ! test -f $SETTINGSFILE
then
    echo "Settings file not found: "$SETTINGSFILE
    exit
fi

if ! test -f $TESTSET
then
    echo "Testset file not found: "$TESTSET
    exit
fi

if ! test -f $BINFILE
then
    echo "SoPlex executable not found: "$BINFILE
    exit
fi

date >$OUTFILE
date >$ERRFILE

# Avoid problems with foreign locales (two separate commands for SunOS)
LANG=C
export LANG

# Determine awk program to use.
AWK=awk
OSTYPE=`uname -s | tr '[:upper:]' '[:lower:]' | sed -e s/cygwin.*/cygwin/ -e s/irix../irix/`

case $OSTYPE in
    osf1)  AWK=gawk ;;
    sunos)  AWK=gawk ;;
    aix)  AWK=gawk ;;
esac

# Create testset
$BINFILE --loadset=$SETTINGSFILE -t$TIME --saveset=$SETFILE

# Solve the instances of the testset
for instance in `cat $TESTSET`
do
    echo @01 $instance
    echo @01 $instance >>$ERRFILE
    $BINFILE $opt --loadset=$SETFILE -c -q -t$TIME $instance 2>>$ERRFILE
    echo =ready=
done | tee -a $OUTFILE
date >>$OUTFILE
date >>$ERRFILE
evaluation.py $SOLUNAME $TESTSET $OUTFILE | tee $RESFILE