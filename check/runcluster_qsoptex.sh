#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#*                                                                           *#
#*                  This file is part of the class library                   *#
#*       SoPlex --- the Sequential object-oriented simPlex.                  *#
#*                                                                           *#
#*  Copyright 1996-2022 Zuse Institute Berlin                                *#
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

# Execute binary and write logfiles, gets called by check_cluster_sqoptex.sh.

# check if tmp-path exists
if test ! -d "${CLIENTTMPDIR}"
then
    echo "Skipping test since the path for the tmp-dir does not exist."
    exit
fi

OUTFILE="${CLIENTTMPDIR}/${BASENAME}.qsoptex.out"
BASFILE="${SOLVERPATH}/results/${BASENAME}.bas"
QSOBASFILE="${SOLVERPATH}/results/${BASENAME}.qsoptex.bas"

# check if basis file exists
if test ! -e "${BASFILE}"
then
    echo "Skipping test since the basis file does not exist."
    exit
fi

uname -a                               > "${OUTFILE}"
echo "@21 ${FILENAME}"                >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date +"@23 %s"                        >> "${OUTFILE}"
sed '/s/^ XU/ XL/g' "${BASFILE}" > "${QSOBASFILE}"
"${EXECNAME}" -B "${QSOBASFILE}" "${INSTANCE}" >> "${OUTFILE}" 2>> "${OUTFILE}"
date +"@24 %s"                        >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"

mv "${OUTFILE}" "${SOLVERPATH}/results/${BASENAME}.qsoptex.out"
