# $Id: check.sh,v 1.2 2001/11/06 23:30:58 bzfkocht Exp $
BINNAME=`basename $2`
OUTFILE=check.$BINNAME.out
ERRFILE=check.$BINNAME.err
RESFILE=check.$BINNAME.result
date >$OUTFILE
date >$ERRFILE
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
        ../$2 $opt $i 2>>$ERRFILE
        echo =ready=
    done
done | tee -a $OUTFILE
date >>$OUTFILE
date >>$ERRFILE
gawk -f check.awk $OUTFILE | tee $RESFILE
 
