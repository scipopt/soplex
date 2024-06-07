#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#*                                                                           *#
#*                  This file is part of the class library                   *#
#*       SoPlex --- the Sequential object-oriented simPlex.                  *#
#*                                                                           *#
#*  Copyright (c) 1996-2024 Zuse Institute Berlin (ZIB)                      *#
#*                                                                           *#
#*  Licensed under the Apache License, Version 2.0 (the "License");          *#
#*  you may not use this file except in compliance with the License.         *#
#*  You may obtain a copy of the License at                                  *#
#*                                                                           *#
#*      http://www.apache.org/licenses/LICENSE-2.0                           *#
#*                                                                           *#
#*  Unless required by applicable law or agreed to in writing, software      *#
#*  distributed under the License is distributed on an "AS IS" BASIS,        *#
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *#
#*  See the License for the specific language governing permissions and      *#
#*  limitations under the License.                                           *#
#*                                                                           *#
#*  You should have received a copy of the Apache-2.0 license                *#
#*  along with SoPlex; see the file LICENSE. If not email soplex@zib.de.     *#
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
