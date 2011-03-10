BINNAME=`basename $2`
TSTNAME=`basename $1 .test`
OUTFILE=cover.$TSTNAME.$BINNAME.out
ERRFILE=cover.$TSTNAME.$BINNAME.err
RESFILE=cover.$TSTNAME.$BINNAME.res
date >$OUTFILE
date >$ERRFILE
for i in `cat $1`
do
    for algo in "" "-e" "-r" "-i" "-e -r" "-e -i" "-i -r" "-e -i -r" 
    do
	for starter in "-c0" "-c1" "-c2"  
	do
	    for simpl in "-s0" "-s1" "-s2" "-s3" "-s4" "-s5"
	    do
		for ratio in "-t0" "-t1" "-t2" 
		do
		    for pricer in "-p0" "-p1" "-p2" "-p3" "-p4" "-p5"
		    do
			opt="-l10 $algo $starter $simpl $ratio $pricer" 

			echo @01 $i ===========
			echo @01 $i =========== >>$ERRFILE
			echo =start= $opt 
			../$2 $opt -l60 $i 2>>$ERRFILE
			echo =ready=
		    done
		done
	    done
	done
    done
done >>$OUTFILE
date >>$OUTFILE
date >>$ERRFILE
gawk -f check.awk $OUTFILE | tee $RESFILE
