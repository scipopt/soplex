# $Id: check.awk,v 1.3 2001/11/08 08:47:18 bzfkocht Exp $
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
    sol["25fv45"] =     5.5018458883E+03;
    sol["80bau3b"] =    9.8723216072E+05;
    sol["adlittle"] =   2.2549496316E+05;
    sol["afiro"] =     -4.6475314286E+02;
    sol["agg"] =       -3.5991767287E+07;
    sol["agg2"] =      -2.0239252356E+07;
    sol["agg3"] =       1.0312115935E+07;
    sol["bandm"] =     -1.5862801845E+02;
    sol["beaconfd"] =   3.3592485807E+04;
    sol["blend"] =     -3.0812149846E+01;
    sol["bnl1"] =       1.9776292856E+03;
    sol["bnl2"] =       1.8112365404E+03;
    sol["boeing1"] =   -3.3521356751E+02;
    sol["boeing2"] =   -3.1501872802E+02;
    sol["bore3d"] =     1.3730803942E+03;
    sol["brandy"] =     1.5185098965E+03;
    sol["capri"] =      2.6900129138E+03;
    sol["cycle"] =     -5.2263930249E+00;
    sol["czprob"] =     2.1851966989E+06;
    sol["d2q06c"] =     1.2278423615E+05;
    sol["d6cube"] =     3.1549166667E+02;
    sol["degen2"] =    -1.4351780000E+03;
    sol["degen3"] =    -9.8729400000E+02;
    sol["dfl001"] =     1.1266396047E+07;
    sol["e226"] =      -1.8751929066E+01;
    sol["etamacro"] =  -7.5571521774E+02;
    sol["fffff800"] =   5.5567961165E+05;
    sol["finnis"] =     1.7279096547E+05;
    sol["fit1d"] =     -9.1463780924E+03;
    sol["fit1p"] =      9.1463780924E+03;
    sol["fit2d"] =     -6.8464293294E+04;
    sol["fit2p"] =      6.8464293232E+04;
    sol["forplan"] =   -6.6421873953E+02;
    sol["ganges"] =    -1.0958636356E+05;
    sol["gfrd-pnc"] =   6.9022359995E+06;
    sol["greenbea"] =  -7.2555248130e+07;
    sol["greenbeb"] =  -4.3022602612E+06;
    sol["grow15"] =    -1.0687094129E+08;
    sol["grow22"] =    -1.6083433648E+08;
    sol["grow7"] =     -4.7787811815E+07;
    sol["israel"] =    -8.9664482186E+05;
    sol["kb2"] =       -1.7499001299E+03;
    sol["lotfi"] =     -2.5264706062E+01;
    sol["maros"] =     -5.8063743701E+04;
    sol["maros-r7"] =   1.4971851665E+06;
    sol["modszk1"] =    3.2061972906E+02;
    sol["nesm"] =       1.4076073035E+07;
    sol["perold"] =    -9.3807580773E+03;
    sol["pilot"] =     -5.5748970615E+02;
    sol["pilot-ja"] =  -6.1131344111E+03;
    sol["pilot-we"] =  -2.7201027439E+06;
    sol["pilot4"] =    -2.5811392641E+03;
    sol["pilot87"] =    3.0171034734E+02
    sol["pilotnov"] =  -4.4972761882E+03;
    sol["qap8"] =       2.0350000000E+02;
    sol["qap12"] =      5.2289435056E+02;
    sol["qap15"] =      1.0409940410E+03;
    sol["recipe"] =    -2.6661600000E+02;
    sol["sc105"] =     -5.2202061212E+01;
    sol["sc205"] =     -5.2202061212E+01;
    sol["sc50a"] =     -6.4575077059E+01;
    sol["sc50b"] =     -7.0000000000E+01;
    sol["scagr25"] =   -1.4753433061E+07;
    sol["scagr7"] =    -2.3313892548E+06;
    sol["scfxm1"] =     1.8416759028E+04;
    sol["scfxm2"] =     3.6660261565E+04;
    sol["scfxm3"] =     5.4901254550E+04;
    sol["scorpion"] =   1.8781248227E+03;
    sol["scrs8"] =      9.0429998619E+02;
    sol["scsd1"] =      8.6666666743E+00;
    sol["scsd6"] =      5.0500000078E+01;
    sol["scsd8"] =      9.0499999993E+02;
    sol["sctap1"] =     1.4122500000E+03;
    sol["sctap2"] =     1.7248071429E+03;
    sol["sctap3"] =     1.4240000000E+03;
    sol["seba"] =       1.5711600000E+04;
    sol["share1b"] =   -7.6589318579E+04;
    sol["share2b"] =   -4.1573224074E+02;
    sol["shell"] =      1.2088253460E+09;
    sol["ship04l"] =    1.7933245380E+06;
    sol["ship04s"] =    1.7987147004E+06;
    sol["ship08l"] =    1.9090552114E+06;
    sol["ship08s"] =    1.9200982105E+06;
    sol["ship12l"] =    1.4701879193E+06;
    sol["ship12s"] =    1.4892361344E+06;
    sol["sierra"] =     1.5394362184E+07;
    sol["stair"] =     -2.5126695119E+02;
    sol["standata"] =   1.2576995000E+03;
    sol["standmps"] =   1.4060175000E+03;
    sol["stocfor1"] =  -4.1131976219E+04;
    sol["stocfor2"] =  -3.9024408538E+04;
    sol["stocfor3"] =  -3.9976661576E+04;
    sol["truss"] =      4.5881584719E+05;
    sol["tuff"] =       2.9214776509E-01;
    sol["vtp"] =        1.2983146246E+05;
    sol["wood1p"] =     1.4429024116E+00;
    sol["woodw"] =      1.3044763331E+00;
}
/=type=/      { type = $2; }
/loading/     { file = $4; }
/rows/        { rows = $3; }
/columns/     { cols = $1; } 
/time/        { time = $4; } 
/iterations/  { iter = $3; }
/value/       { obj  = $4; }
/=start=/     {
   type = "";
   for(i = 2; i <= NF; i++)
      type = type substr($i, 2);
}
/ready/       {
    n = split(file, a, "/");
    split(a[n], b, ".");
    name = b[1];

    abserr = abs(sol[name] - obj);

    if (sol[name] == "")
	print name, "nicht gefunden";
    else
       relerr = abserr / abs(sol[name]) * 1000.0

    if (name == prevname)
        printf("%23s", "");
    else
    {
        printf("----------------------------------------------------------------------------\n");
	printf("%-10s %5d %5d ", name, rows, cols);
    }
    printf("%s %5d %8.2f %+e ", type, iter, time, obj);

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

    sum[type] += time;
    cnt[type]++;
    counts++;
    times += time;
    
    prevname = name;
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



