/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: lpconv.cpp,v 1.1 2002/03/10 10:00:59 bzfkocht Exp $"

#include <assert.h>
#include <iostream>
#include <fstream>

#include "spxdefines.h"
#include "spxlp.h"

using namespace soplex;

int main(int argc, char **argv)
{
   const char* banner =
   "************************************************************************\n"
   "*                                                                      *\n"
   "*       LPConv --- Convert LPF to MPS format.                          *\n"
   "*                  Release 1.0.0                                       *\n"
   "*    Copyright (C) 2002 Konrad-Zuse-Zentrum                            *\n"
   "*                       fuer Informationstechnik Berlin                *\n"
   "*                                                                      *\n"
   "*  LPConv is distributed under the terms of the ZIB Academic Licence.  *\n"
   "*  You should have received a copy of the ZIB Academic License         *\n"
   "*  along with SoPlex; If not email to soplex@zib.de.                   *\n"
   "*                                                                      *\n"
   "************************************************************************\n"
   ;

   const char* usage =
   "[options] input-file output-file\n\n"
   "          input-file can be either in MPS or LPF format\n\n"
   "options:  (*) indicates default\n" 
   " -v        show program version\n"
   " -h        show this help\n"
   ;

   int verbose = 1;

   for(optidx = 1; optidx < argc; optidx++)
   {
      if (*argv[optidx] != '-')
         break;

      switch(argv[optidx][1])
      {
      case 'v' :
         std::cout << banner << std::endl;
         exit(0);
      case 'h' :
      case '?' :
         std::cout << banner << std::endl;
         /*FALLTHROUGH*/
      default :
         std::cerr << "usage: " << argv[0] << " " << usage << std::endl;
         exit(0);
      }
   }
   if ((argc - optidx) < 2)
   {
      std::cerr << argv[0] << ":" << usage << std::endl;
      exit(0);
   }
   const char* inpfile  = argv[optidx];
   const char* outfile  = argv[optidx + 1];

   Param::setVerbose(verbose);

   SPxLP       lp;
   NameSet     rownames;
   NameSet     colnames;
   DIdxSet     intvars;

   std::ifstream ifile(inpfile);

   if (!ifile)
   {
      std::cerr << "Can't open file: " << inpfile << std::endl;
      exit(1);
   }
   if (!lp.read(ifile, rownames, colnames, intvars))
   {
      std::cerr << "Error while reading file: " << inpfile << std::endl;
      exit(1);
   }

   std::ofstream ofile(outfile);
   
   if (!ofile)
   {
      std::cerr << "Can't open file: " << outfile << std::endl;
      exit(1);
   }
   lp.writeMPS(ofile, rownames, colnames, intvars);

   return 0;
}




