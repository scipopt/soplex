#!/usr/bin/env bash

# solve a given testset with given settings and time limit
# parameters:
# 1: name of testset (has to be in check/testset)
# 2: path to soplex executable
# 3: name of settings (has to be in settings)
# 4: time limit
# 5: results directory

TEST=$1
TESTSET=testset/$TEST.test
SOLUNAME=testset/$TEST.solu

BINFILE=$2
BINNAME=`basename $2`

# get host name
HOST=`uname -n | sed 's/\(.zib.de\)//g'`
BINNAME=$BINNAME.$HOST

SETTINGS=$3
SETTINGSFILE=../settings/$SETTINGS.set

TIME=$4

RESDIR=$5

OUTFILE=$RESDIR/check.$TEST.$BINNAME.$SETTINGS.out
ERRFILE=$RESDIR/check.$TEST.$BINNAME.$SETTINGS.err
RESFILE=$RESDIR/check.$TEST.$BINNAME.$SETTINGS.res
SETFILE=$RESDIR/check.$TEST.$BINNAME.$SETTINGS.set

# create results directory
mkdir -p $RESDIR

# Abort if files are missing
if ! test -f $SETTINGSFILE
then
    if [ "$SETTINGSFILE" == "../settings/default.set" ];
    then
        touch $SETTINGSFILE
    else
        echo "Settings file not found: "$SETTINGSFILE
        exit 1
    fi
fi

if ! test -f $TESTSET
then
    echo "Testset file not found: "$TESTSET
    exit 1
fi

if ! test -f $BINFILE
then
    echo "SoPlex executable not found: "$BINFILE
    exit 1
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
    $BINFILE $opt --loadset=$SETFILE -v4 --int:displayfreq=10000 -c -q -t$TIME $instance 2>>$ERRFILE
    echo =ready=
done | tee -a $OUTFILE
date >>$OUTFILE
date >>$ERRFILE

# check whether python is available
if command -v python >/dev/null 2>&1
then
	python evaluation.py $OUTFILE | tee $RESFILE
else
	./evaluation.sh $OUTFILE | tee $RESFILE
fi
