#!/bin/sh
PATH=$PATH:/client/bin
export PATH
cd /zibis/Optimization/Software/Soplex
metasend -z -b -S 1000000 -t $1 -F koch@zib.de -s Soplex-1.2.1 -e quoted-printable -m text/plain -f greeting.txt -D Instructions -n -e base64 -f soplex-1.2.1.tar.gz -m application/x-gzip -D soplex-1.2.1.tar.gz 

