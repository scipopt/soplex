# $Id: check.sh,v 1.6 2001/11/22 16:30:00 bzfkocht Exp $
BINNAME=`basename $2`
TSTNAME=`basename $1 .test`
OUTFILE=check.$TSTNAME.$BINNAME.out
ERRFILE=check.$TSTNAME.$BINNAME.err
RESFILE=check.$TSTNAME.$BINNAME.res
date >$OUTFILE
date >$ERRFILE
case $TSTNAME in
mittelbach) timelimit="-l10000" ;;
*)          timelimit="" ;;
esac
for i in `cat $1`
do
    echo @01 $i ===========
    echo @01 $i =========== >>$ERRFILE
    for k in 1 2 3 4
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
        esac
        ../$2 $opt -s3 $timelimit $i 2>>$ERRFILE
        echo =ready=
    done
done | tee -a $OUTFILE
date >>$OUTFILE
date >>$ERRFILE
gawk -f check.awk $OUTFILE | tee $RESFILE
 
