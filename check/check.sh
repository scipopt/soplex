# $Id: check.sh,v 1.19 2004/03/22 11:35:23 bzfpfend Exp $
# Parameters
# $1 Name of the test, e.g. netlib (needs netlib.test, netlib.solu)
# $2 Path/Name of the binary, e.g. ../bin/soplex.linux.x86.gnu.opt
# $3 Algorithms to test (1...6), e.g. "1 2 3 4"
# $4 Limits, e.g. -l10000 as time limit.
BINNAME=`basename $2`
TSTNAME=`basename $1 .test`
OUTFILE=check.$TSTNAME.$BINNAME.out
ERRFILE=check.$TSTNAME.$BINNAME.err
RESFILE=check.$TSTNAME.$BINNAME.res
date >$OUTFILE
date >$ERRFILE
#
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
        esac
        $2 $opt -q -s0 $4 $i 2>>$ERRFILE
        echo =ready=
    done
done | tee -a $OUTFILE
date >>$OUTFILE
date >>$ERRFILE
gawk -f check.awk $TSTNAME.solu $OUTFILE | tee $RESFILE
 
