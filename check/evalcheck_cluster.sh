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

# evaluate logfiles from testrun started by 'make testcluster'
#
# USAGE, from directory check, call
#              ./evalcheck_cluster.sh results/check.*.eval

export LANG=C

REMOVE=0
DECOMP=0
AWKARGS=""
FILES=""

for i in $@
do
    if test ! -e "${i}"
    then
        if test "${i}" = "-r"
        then
            REMOVE=1
        else
            AWKARGS="${AWKARGS} ${i}"
        fi

        if test "${i}" = "-d"
        then
            DECOMP=1
        fi
    else
        FILES="${FILES} ${i}"
    fi
done

for FILE in ${FILES}
do
    DIR=$(dirname "${FILE}")
    EVALFILE=$(basename "${FILE}" .eval)
    EVALFILE=$(basename "${EVALFILE}" .out)

    OUTFILE="${DIR}/${EVALFILE}.out"
    ERRFILE="${DIR}/${EVALFILE}.err"
    RESFILE="${DIR}/${EVALFILE}.res"
    SETFILE="${DIR}/${EVALFILE}.set"

    # check if the eval file exists; if this is the case construct the overall solution files
    if test -e "${DIR}/${EVALFILE}.eval"
    then
        # in case an output file exists, copy it away to save the results
        DATEINT=$(date +"%s")
        if test -e "${OUTFILE}"
        then
            cp "${OUTFILE}" "${OUTFILE}.old-${DATEINT}"
        fi
        if test -e "${ERRFILE}"
        then
            cp "${ERRFILE}" "${ERRFILE}.old-${DATEINT}"
        fi

        echo > "${OUTFILE}"
        echo > "${ERRFILE}"
        echo "create overall output and error file for ${EVALFILE}"

        for i in $(cat "${DIR}/${EVALFILE}.eval") DONE
        do
            if test "${i}" = "DONE"
            then
                break
            fi

            FILE="${i}.out"
            if test -e "${FILE}"
            then
                cat "${FILE}" >> "${OUTFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            else
                echo "Missing ${i}"
            fi

            FILE="${i}.perplex.out"
            if test -e "${FILE}"
            then
                echo             >> "${OUTFILE}"
                echo "=perplex=" >> "${OUTFILE}"
                cat "${FILE}" >> "${OUTFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            fi

            FILE="${i}.qsoptex.out"
            if test -e "${FILE}"
            then
                echo             >> "${OUTFILE}"
                echo "=qsoptex=" >> "${OUTFILE}"
                cat "${FILE}" >> "${OUTFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            fi

            echo           >> "${OUTFILE}"
            echo "=ready=" >> "${OUTFILE}"

            FILE="${i}.err"
            if test -e "${FILE}"
            then
                cat "${FILE}" >> "${ERRFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            fi

            FILE="${i}.set"
            if test -e "${FILE}"
            then
                cp "${FILE}" "${SETFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            fi

            if test "${REMOVE}" = "1"
            then
                rm -f "${i}.bas"
            fi
        done

        if test "${REMOVE}" = "1"
        then
            rm -f "${DIR}/${EVALFILE}.eval"
        fi
    fi

    # check if the out file exists
    if test -e "${DIR}/${EVALFILE}.out"
    then
        echo "create results for ${EVALFILE}"

        if test "${DECOMP}" = "1"
        then
            ./evaluation.py -d "${OUTFILE}" | tee "${RESFILE}"
        else
            ./evaluation.py "${OUTFILE}" | tee "${RESFILE}"
        fi
    fi
done
