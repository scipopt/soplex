# $Id: check.sh,v 1.12 2002/02/14 15:07:42 bzfkocht Exp $
BINNAME=`basename $2`
TSTNAME=`basename $1 .test`
OUTFILE=check.$TSTNAME.$BINNAME.out
ERRFILE=check.$TSTNAME.$BINNAME.err
RESFILE=check.$TSTNAME.$BINNAME.res
date >$OUTFILE
date >$ERRFILE
case $TSTNAME in
mittelmann|zib) timelimit="-l10000" 
            algorithm="1 2" ;;
*)          timelimit="" 
            algorithm="1 2 3 4 5 6" ;;
esac
for i in `cat $1`
do
    echo @01 $i ===========
    echo @01 $i =========== >>$ERRFILE
    for k in $algorithm
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
        ../$2 $opt $timelimit $i 2>>$ERRFILE
        echo =ready=
    done
done | tee -a $OUTFILE
date >>$OUTFILE
date >>$ERRFILE
gawk -f check.awk $TSTNAME.solu $OUTFILE | tee $RESFILE
 
