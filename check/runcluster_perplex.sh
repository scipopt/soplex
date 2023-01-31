#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#*                                                                           *#
#*                  This file is part of the class library                   *#
#*       SoPlex --- the Sequential object-oriented simPlex.                  *#
#*                                                                           *#
#*  Copyright (c) 1996-2023 Zuse Institute Berlin (ZIB)                      *#
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

# Execute binary and write logfiles, gets called by check_cluster_perplex.sh.

# check if tmp-path exists
if test ! -d "${CLIENTTMPDIR}"
then
    echo Skipping test since the path for the tmp-dir does not exist.
    exit
fi

OUTFILE="${CLIENTTMPDIR}/${BASENAME}.perplex.out"
BASFILE="${SOLVERPATH}/results/${BASENAME}.bas"

# check if basis file exists
if test ! -e "${BASFILE}"
then
    echo Skipping test since the basis file does not exist.
    exit
fi

uname -a                              > "${OUTFILE}"
echo "@11 ${FILENAME}"                >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date +"@13 %s"                        >> "${OUTFILE}"
"${EXECNAME}" "${INSTANCE}" "${BASFILE}" >> "${OUTFILE}" 2>> "${OUTFILE}"
date +"@14 %s"                        >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"

mv "${OUTFILE}" "${SOLVERPATH}/results/${BASENAME}.perplex.out"
