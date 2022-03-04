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
