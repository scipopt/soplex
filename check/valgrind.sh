#! /bin/bash
#
# Ths script runs the binary under control of valgrind.
#
# Parameters
# $1 Name of the test, e.g. netlib (needs netlib.test, netlib.solu)
# $2 Path/Name of the binary, e.g. ../bin/soplex.linux.x86.gnu.opt
# $3 Algorithms to test (1...6), e.g. "1 2 3 4"
# $4 Limits, e.g. -l10000 as time limit.
# $5 valgrind command with parameters (inclusive valgrind binary name)
# $6 potential suppression file name
BINNAME=`basename $2`
TSTNAME=`basename $1 .test`
OUTFILE=check.$TSTNAME.$BINNAME.out
ERRFILE=check.$TSTNAME.$BINNAME.err
RESFILE=check.$TSTNAME.$BINNAME.res
date >$OUTFILE
date >$ERRFILE

# echo '$1->' $1
# echo '$2->' $2
# echo '$3->' $3
# echo '$4->' $4
# echo '$5->' $5
# echo '$6->' $6



VALGRINDCMD=$5
VSUPPNAME=$6

# look for an error suppression file
if [ -f $VSUPPNAME ]
then
  VALGRINDCMD="$VALGRINDCMD --suppressions=$VSUPPNAME"
fi

for i in `cat $1`
do
    echo @01 $i ===========
    echo @01 $i =========== >>$ERRFILE
    for k in $3
    do
        case $k in
	1)  echo =type= LC
	    opt="" ;;
	2)  echo =type= EC
	    opt="-e" ;;
	3)  echo =type= LR
	    opt="-r" ;;
	4)  echo =type= ER
            opt="-e -r" ;;
	5)  echo =type= LCi
	    opt="-i" ;;
	6)  echo =type= ECi
            opt="-e -i" ;;
	7)  echo =type= LCd
	    opt="-p2" ;;
	8)  echo =type= ECd
	    opt="-e -p2" ;;
	9)  echo =type= LCh
	    opt="-t1" ;;
	10) echo =type= ECh
	    opt="-e -t1" ;;
	11) echo =type= LCm
	    opt="-p1" ;;
	12) echo =type= ECm
	    opt="-e -p1" ;;
        esac

        $VALGRINDCMD $2 $opt -q $4 $i 2>>$ERRFILE
        echo =ready=
    done
done | tee -a $OUTFILE
date >>$OUTFILE
date >>$ERRFILE
gawk -f check.awk $TSTNAME.solu $OUTFILE | tee $RESFILE
 
