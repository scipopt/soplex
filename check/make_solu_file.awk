# $Id: make_solu_file.awk,v 1.4 2006/02/27 17:31:51 bzfhille Exp $
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
      printf( "=opt= %-16s %-40s\n", file, sol_cplex[ file ] );

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






