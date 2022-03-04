#!/bin/bash
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

# Call with 'make test' and 'make check' (via check.sh)
#

# solve a given testset with given settings and time limit
# parameters:

TSTNAME="${1}"    # name of testset (has to be in check/testset)
EXECUTABLE="${2}" # path to soplex executable
SETTINGS="${3}"   # name of settings (has to be in settings)
TIME="${4}"       # time limit
OUTPUTDIR="${5}"  # results directory


# check if all variables defined (by checking the last one)
if test -z "${OUTPUTDIR}"
then
    echo Skipping test since not all variables are defined
    echo "TSTNAME       = ${TSTNAME}"
    echo "EXECUTABLE    = ${EXECUTABLE}"
    echo "SETTINGS      = ${SETTINGS}"
    echo "TIME          = ${TIME}"
    echo "OUTPUTDIR     = ${OUTPUTDIR}"
    exit 1;
fi

# call routines for creating the result directory, checking for existence
# of passed settings, etc
# defines the following environment variables: SOPLEXPATH, FULLTSTNAME, SOLUFILE, SETTINGSFILE
. ./configuration_set.sh

BINNAME=$(basename "${EXECUTABLE}")
# get host name
HOST=$(uname -n | sed 's/\(.zib.de\)//g')
BINID="${BINNAME}.${HOST}"

OUTFILE="${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${SETTINGS}.out"
ERRFILE="${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${SETTINGS}.err"
RESFILE="${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${SETTINGS}.res"
SETFILE="${OUTPUTDIR}/check.${TSTNAME}.${BINID}.${SETTINGS}.set"

# create results directory
mkdir -p "${OUTPUTDIR}"

if ! test -f "${EXECUTABLE}"
then
    echo "SoPlex executable not found: ${EXECUTABLE}"
    exit 1
fi

date >"${OUTFILE}"
date >"${ERRFILE}"

# Avoid problems with foreign locales (two separate commands for SunOS)
LANG=C
export LANG

# Determine awk program to use.
AWK=awk
OSTYPE=$(uname -s | tr '[:upper:]' '[:lower:]' | sed -e s/cygwin.*/cygwin/ -e s/irix../irix/)

case "${OSTYPE}" in
    osf1)  AWK=gawk ;;
    sunos)  AWK=gawk ;;
    aix)  AWK=gawk ;;
esac

# Create testset
"${EXECUTABLE}" --loadset="${SETTINGSFILE}" -t"${TIME}" --saveset="${SETFILE}"

# Solve the instances of the testset
for instance in $(cat "${FULLTSTNAME}")
do
    echo "@01 ${instance}"
    echo "@01 ${instance}" >> "${ERRFILE}"
    ${EXECUTABLE} --loadset="${SETFILE}" -v4 --int:displayfreq=10000 -c -q -t"${TIME}" "${instance}" 2>> "${ERRFILE}"
    echo "=ready="
done | tee -a "${OUTFILE}"
date >> "${OUTFILE}"
date >> "${ERRFILE}"

# check whether python is available
if command -v python >/dev/null 2>&1
then
    python evaluation.py "${OUTFILE}" | tee "${RESFILE}"
else
    ./evaluation.sh "${OUTFILE}" | tee "${RESFILE}"
fi
