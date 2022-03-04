#!/bin/sh
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

#
# Shell script for generating .solu file from a .test file using CPLEX and perplex.
#
# .test file may be generated by editing the output of the following command:
# ls -1 --color=never --indicator-style=none -d LP/kennington/*
#

if [ $# -ne 1 ]
then
    echo "usage: ${0} <instance-file>"
    exit 0
fi


INSTANCE_LIST="${1}"
PERPLEX_LOG="${INSTANCE_LIST}.perplex"
CPLEX_LOG="${INSTANCE_LIST}.cplex"
BASIS_DIR=basis
ERROR=0

# Checks for sanity.
if (test -z $(which cplex 2> /dev/null))
then
    echo "CPLEX not found"
    ERROR=1
fi

if (test -z $(which perplex 2> /dev/null))
then
    echo "perplex not found"
    ERROR=1
fi

if !(test -e "${INSTANCE_LIST}")
then
    echo "file ${INSTANCE_LIST} not found"
    ERROR=1
fi

if (test "${ERROR}" -gt 0)
then
    exit 1
fi


#
# Main part.
#

rm -f "${PERPLEX_LOG}" "${CPLEX_LOG}"
mkdir -p basis

# CPLEX and perplex loop
for file in $(cat "${INSTANCE_LIST}");
do
    # Setup CPLEX batch file depending on file format.
    format_known=1

    case "${file}" in
        *.mps)    # MPS format
            ext="mps"
            ;;
        *.mps.gz) # gzipped MPS format
            ext="mps.gz"
            ;;
        *.lp)     # LP format
            ext="lp"
            ;;
        *.lp.gz)  # gzipped LP format
            ext="lp.gz"
            ;;
        *)        # unknown format
            format_known=0
            ;;
    esac

    if [ "${format_known}" -gt 0 ]; then
        base=$(basename "${file}" ".${ext}")

        # Solve using CPLEX.
        cplex <<-EOF
    set logfile ${CPLEX_LOG}
    set simplex tol feas 1e-9
    set simplex tol opt 1e-9
    r ${file}
    change problem lp
    opt
    w work.bas
    quit
EOF

        # Name the .bas file correctly, run perplex and write everything to
        # logfile for later evaluation.
        echo "@01 ${file} ===========" >> "${PERPLEX_LOG}"

        if [ -e work.bas ]; then
            echo "work.bas -> ${base}.bas.gz"
            mv -f work.bas "${BASIS_DIR}/${base}.bas"
            rm -f "${BASIS_DIR}/${base}.bas.gz"
            gzip "${BASIS_DIR}/${base}.bas"

            if ( test "${ext}" == "lp" ) || ( test "${ext}" == "lp.gz" )
            then
                cplex <<-EOF
    r "${file}"
    w "${base}.mps.gz"
    quit
EOF
                file="${base}.mps.gz"
            fi

            perplex "${file}" "${BASIS_DIR}/${base}.bas.gz" | tee -a "${PERPLEX_LOG}"

            if ( test "${ext}" == "lp" ) || ( test "${ext}" == "lp.gz" )
            then
                rm -f "${file}"
            fi
        else
            echo "some error in CPLEX run" >> "${PERPLEX_LOG}"
        fi
    fi
done

# Evaluate perplex log via awk script.
base=$(basename "${INSTANCE_LIST}" .test)
awk -f make_solu_file.awk "${INSTANCE_LIST}.cplex" "${INSTANCE_LIST}.perplex" | tee "${base}.solu"

exit
