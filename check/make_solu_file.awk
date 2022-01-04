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
# Takes two files:  <CPLEX log> <Perplex log> (order is critical!)
#
BEGIN {
    first = 1;
    feasible = 0;
    optimal = 0;
    obj = "";
    file = "";
}

# Read CPLEX log and associate objective value to instance name.
/Problem '/ {
    file_cplex = $2;
    gsub( "'", "", file_cplex );
    sub( /.*\//, "", file_cplex );

    # Eliminate suffix if present.
    sub( /\.mps/, "", file_cplex );
    sub( /\.lp/, "", file_cplex );
    sub( /\.gz/, "", file_cplex );
}

# Must be before /Objective/ to ensure correct infeasibility detection (CPLEX may
# give an objective value although it already found infeasibility).
/Infeasible/ {
    obj_cplex = "infeasible";
}

/Objective =/ {
    obj_cplex = $7;
}

/Basis written/ {
    #  printf( "%s: %s\n", file_cplex, obj_cplex );

    sol_cplex[ file_cplex ] = obj_cplex;
    file_cplex = "";
    obj_cplex = "E";
}

# Read Perplex log.
/@01/ {
    if ( !first )
        if ( feasible && optimal )
        printf( "=opt= %-16s %-40s\n", file, obj );
    else
        printf( "=opt= %-16s %-40s # CPLEX not optimal/feasible\n", file, sol_cplex[ file ] );

    first = 0;
    feasible = 0;
    optimal = 0;
    obj = "";
    file = "";
}
/Name:/ {
    file = $2;
    # Eliminate suffix if present.
    sub( /\.mps/, "", file );
    sub( /\.lp/, "", file );
    sub( /\.gz/, "", file );
}
/Solution is optimal/          { optimal = 1; }
/Solution is feasible/         { feasible = 1; }
/objective/                    { obj = $3; }


END {
    if ( feasible && optimal )
        printf( "=opt= %-16s %-40s\n", file, obj );
    else
        printf( "=opt= %-16s %-40s\n", file, sol_cplex[ file ] );
}
