# $Id: check.awk,v 1.14 2002/02/20 16:30:35 bzfkocht Exp $
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*   File....: check.awk                                                     *
#*   Name....: SoPlex Check Report Generator                                 *
#*   Author..: Thorsten Koch                                                 *
#*   Copyright by Author, All rights reserved                                *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
function abs(x)
{
    return x < 0 ? -x : x;
}
BEGIN {
    print "$Id: check.awk,v 1.14 2002/02/20 16:30:35 bzfkocht Exp $";
    print "";
}
/=opt=/          { sol[$2] = $3; }
/=type=/         { type = $2; }
/loading/        { file = $4; }
/rows/           { rows = $3; }
/columns/        { cols = $1; } 
/solution time/  { time = $4; } 
/iterations/     { iter = $3; }
/value/          { obj  = $4; }
/infeasible/     { infeas = 1; } 
/unbounded/      { infeas = 1; }
/time limit/     { timeout = 1; }
/=start=/        {
   type = "";
   for(i = 2; i <= NF; i++)
      type = type substr($i, 2);
}
/ready/       {
    n = split(file, a, "/");
    split(a[n], b, ".");
    name = b[1];

    if (sol[name] == "")
       print name, "nicht gefunden";
    else
    {
        if (name == prevname)
            printf("%25s", "");
        else
        {
            printf("------------------------------------------------------------------------------\n");
	    printf("%-10s %6d %6d ", name, rows, cols);
        }
	printf("%-3s %7d %8.2f ", type, iter, time);

        if (infeas)
	    printf("%-14s", "infeasible");
	else if (timeout)
	    printf("%-14s", "timeout");
	else
	    printf("%+e ", obj);

	if (timeout)
	   printf("\n");
	else
	{
            if (!infeas && sol[name] != "infeasible")
            {
                abserr = abs(sol[name] - obj);

                if (sol[name] != 0.0)
                    relerr = abserr / abs(sol[name]) * 1000.0
                else
                    relerr = 0.0;

    	        if ((abserr < 1e-6) || (relerr < 0.01))
		{
		   printf("ok\n");
		   pass[type]++;
		   passes++;
		}
		else
		{
		   printf("error %e\n", abserr);
		   fail[type]++;
		   fails++;
		}
	    }
	    else
	    {
	       if (infeas == 1 && sol[name] == "infeasible")
	       {
		  printf("ok\n");
		  pass[type]++;
		  passes++;
	       }
	       else
	       {
		  if (infeas && sol[name] != "infeasible")
		     printf("error %e\n", abs(sol[name]));
		  else
		     printf("error infeasible\n");
		  
		  fail[type]++;
		  fails++;
	       }
	    }
	}
        sum[type] += time;
        cnt[type]++;
        counts++;
        times += time;
    }
    prevname = name;
    timeout  = 0;
    infeas   = 0;
    obj      = 0;
    iter     = 0;
    time     = 0;
    rows     = 0;
    cols     = 0;
}
END {
    printf("\n------------------------------------------\n");
    printf("Alg            Cnt  Pass  Fail       Time\n");
    printf("------------------------------------------\n");
    for(i in sum)
    {
	printf("%-12s %5d %5d %5d %10.2f \n", 
	   i, cnt[i], pass[i], fail[i], sum[i]);
    }
    printf("------------------------------------------\n");
    printf("%-12s %5d %5d %5d %10.2f \n",
       "Sum", counts, passes, fails, times);
    printf("------------------------------------------\n");
}



