# $Id: make_solu_file.awk,v 1.1 2005/11/18 09:49:01 bzfhille Exp $
BEGIN {
  first = 1;
  feasible = 0;
  optimal = 0;
  obj = "";
  file = "";
}
/@01/ {
  if ( !first ) 
    if ( feasible && optimal )
      printf( "=opt= %-16s %-40s\n", file, obj );

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
#    printf( "%-12s %1d %1d %-40s\n", file, feasible, optimal, obj );
}






