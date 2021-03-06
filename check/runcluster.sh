#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#*                                                                           *#
#*                  This file is part of the class library                   *#
#*       SoPlex --- the Sequential object-oriented simPlex.                  *#
#*                                                                           *#
#*    Copyright (C) 1996-2022 Konrad-Zuse-Zentrum                            *#
#*                            fuer Informationstechnik Berlin                *#
#*                                                                           *#
#*  SoPlex is distributed under the terms of the ZIB Academic Licence.       *#
#*                                                                           *#
#*  You should have received a copy of the ZIB Academic License              *#
#*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  *#
#*                                                                           *#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#

# Execute binary and write logfiles, gets called by check_cluster.sh.

# check if tmp-path exists
if test ! -d "${CLIENTTMPDIR}"
then
    echo Skipping test since the path for the tmp-dir does not exist.
    exit
fi

OUTFILE="${CLIENTTMPDIR}/${BASENAME}.out"
ERRFILE="${CLIENTTMPDIR}/${BASENAME}.err"
BASFILE="${CLIENTTMPDIR}/${BASENAME}.bas"
TMPFILE="${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.tmp"

uname -a                               > "${OUTFILE}"
uname -a                               > "${ERRFILE}"

# function to copy back the ${OUTPUTDIR} and delete temporary files
function cleanup {
    rm -f "${TMPFILE}"

    mv "${OUTFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.out"
    mv "${ERRFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.err"
    mv "${BASFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.bas"
}

# ensure TMPFILE is deleted and ${OUTPUTDIR} are copied when exiting (normally or due to abort/interrupt)
trap cleanup EXIT

echo "@01 ${FILENAME}"                >> "${OUTFILE}"
echo "@01 ${FILENAME}"                >> "${ERRFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${OUTFILE}"
date                                  >> "${ERRFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date +"@03 %s"                        >> "${OUTFILE}"
"${EXECNAME}" --loadset="${SOLVERPATH}/../settings/${SETTINGS}.set" --writebas="${BASFILE}" -v4 --int:displayfreq=10000 -c -q -t"${TIMELIMIT}" "${INSTANCE}" >> "${OUTFILE}" 2>> "${ERRFILE}"
date +"@04 %s"                        >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${ERRFILE}"

rm -f "${TMPFILE}"
