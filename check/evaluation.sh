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

# Very simple evaluation script,
#
# Called by 'test.sh'

OUTFILE="${1}"

SOPLEXOPTIMAL=$(grep -c "SoPlex status       : problem is solved \[optimal\]" "${OUTFILE}")
optimal=$((SOPLEXOPTIMAL / 2))
total=$(grep -c "=ready=" "${OUTFILE}")

echo
echo 'Summary:'
echo
echo 'total  : '${total}
echo 'optimal: '${optimal}
echo
echo 'please install python for better statistics and evaluation'
