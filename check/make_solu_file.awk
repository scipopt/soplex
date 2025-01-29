#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#*                                                                           *#
#*                  This file is part of the class library                   *#
#*       SoPlex --- the Sequential object-oriented simPlex.                  *#
#*                                                                           *#
#*    Copyright (c) 1996-2025 Zuse Institute Berlin (ZIB)                    *#
#*  Copyright (c) 1996-2025 Zuse Institute Berlin (ZIB)                      *#
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
