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
