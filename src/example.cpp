/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2001 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: example.cpp,v 1.19 2002/01/13 10:12:57 bzfkocht Exp $"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "spxsolver.h"

#include "timer.h"
#include "spxpricer.h"
#include "spxdefaultpr.h"
#include "spxparmultpr.h"
#include "spxdevexpr.h"
#include "spxhybridpr.h"
#include "spxsteeppr.h"
#include "spxweightpr.h"
#include "spxratiotester.h"
#include "spxharrisrt.h"
#include "spxdefaultrt.h"
#include "spxfastrt.h"
#include "spxsimplifier.h"
#include "spxaggregatesm.h"
#include "spxredundantsm.h"
#include "spxrem1sm.h"
#include "spxscale.h"
#include "spxgeneralsm.h"
#include "spxsumst.h"
#include "spxweightst.h"
#include "spxvectorst.h"
#include "slufactor.h"

using namespace soplex;

/** Here comes a simple derived class from #SoPlex, which uses #terminate() as
 *  callback method for outputting statistics.
 */
class MySoPlex : public SPxSolver
{
private:
   SLUFactor m_slu;

public:
   /// default constructor
   MySoPlex(Type p_type = LEAVE, Representation p_rep = COLUMN)
      : SPxSolver(p_type, p_rep)
   {}

   virtual bool terminate()
   {
      if (iteration() % 100 == 0)
         std::cout << iteration() << ":\t" << value() << std::endl;

      return SoPlex::terminate();
   }

   void setUtype(SLUFactor::UpdateType tp)
   {
      m_slu.setUtype(tp);
   }
};

/**@todo Flag to print the solution variables with names. */
int main(int argc, char **argv)
{
   const char* banner =
   "************************************************************************\n"
   "*                                                                      *\n"
   "*       SoPlex --- the Sequential object-oriented simPlex.             *\n"
   "*                  Release 1.2.0                                       *\n"
   "*    Copyright (C) 1997-1999 Roland Wunderling                         *\n"
   "*                  1997-2001 Konrad-Zuse-Zentrum                       *\n"
   "*                            fuer Informationstechnik Berlin           *\n"
   "*                                                                      *\n"
   "*  SoPlex is distributed under the terms of the ZIB Academic Licence.  *\n"
   "*  You should have received a copy of the ZIB Academic License         *\n"
   "*  along with SoPlex; If not email to soplex@zib.de.                   *\n"
   "*                                                                      *\n"
   "************************************************************************\n"
   ;

   const char* usage =
   "[options] file\n\n"
   "          file can be either in MPS or LPF format\n\n"
   "options: (*) indicates default\n"
   " -e       select entering algorithm (default is leaving)\n"
   " -r       select row wise representation (default is column)\n"
   " -i       select Eta-update (default is Forest-Tomlin)\n"
   " -lSec    set timlimit to Sec seconds\n"
   " -dDelta  set maximal allowed bound violation to Delta (1e-6)\n"
   " -zZero   set zero tolerance to Zero (1e-16)\n\n"
   "Simplifier:         Starter:        Pricer:           Ratiotester:\n"
   " -s0  none           -c0  none*      -p0  Textbook     -t0  Textbook\n" 
   " -s1  General        -c1  Weight     -p1  ParMult      -t1  Harris\n"
   " -s2  Aggregate      -c2  Sum        -p2  Devex        -t2  Fast*\n"
   " -s3  Remove-1*      -c3  vector     -p3  Hybrid\n"
   " -s4  Redundant                      -p4  Steep*\n"
   " -s5  Scale                          -p5  Weight\n" 
   ;

   const char*            filename;
   SoPlex::Type           type           = SoPlex::LEAVE;
   SoPlex::Representation representation = SoPlex::COLUMN;
   SLUFactor::UpdateType  update         = SLUFactor::FOREST_TOMLIN;
   SPxSimplifier*         simplifier     = 0;
   SPxStarter*            starter        = 0;
   SPxPricer*             pricer         = 0;
   SPxRatioTester*        ratiotester    = 0;
   NameSet                rownames;
   NameSet                colnames;
   int                    starting       = 0;
   int                    pricing        = 4;
   int                    ratiotest      = 2;
   int                    simplifing     = 3;
   double                 timelimit      = -1.0;
   double                 delta          = 1e-6;
   double                 epsilon        = 1e-16;
   int                    optind;

   for(optind = 1; optind < argc; optind++)
   {
      if (*argv[optind] != '-')
         break;

      switch(argv[optind][1])
      {
      case 'c' :
         starting = atoi(&argv[optind][2]);
         break;
      case 'd' :
         delta = atof(&argv[optind][2]);
         break;
      case 'e':
         type = SoPlex::ENTER;
         break;
      case 'i' :
         update = SLUFactor::ETA;
         break;
      case 'l' :
         timelimit = atof(&argv[optind][2]);
         break;
      case 'p' :
         pricing = atoi(&argv[optind][2]);
         break;
      case 'r' :
         representation = SoPlex::ROW;
         break;
      case 's' :
         simplifing = atoi(&argv[optind][2]);
         break;
      case 't' :
         ratiotest = atoi(&argv[optind][2]);
         break;
      case 'v' :
         std::cout << banner << std::endl;
         exit(0);
      case 'z' :
         epsilon = atof(&argv[optind][2]);
         break;
      case 'h' :
      case '?' :
         std::cerr << "usage: " << argv[0] << " " << usage << std::endl;
         exit(0);
      default :
         abort();
      }
   }
   if ((argc - optind) < 1)
   {
      std::cerr << argv[0] << ":" << usage << std::endl;
      exit(0);
   }
   filename = argv[optind];

   std::cout.setf(std::ios::scientific | std::ios::showpoint);

   MySoPlex work(type, representation);

   work.setUtype(update);
   work.setTermination(timelimit);
   work.setDelta(delta);
   work.setEpsilon(epsilon);

   assert(work.isConsistent());

   std::cout << (type == SoPlex::ENTER ? "Entering" : "Leaving")
             << " algorithm" 
             << std::endl
             << (representation == SoPlex::ROW ? "Row" : "Column")
             << " representation" 
             << std::endl
             << (update == SLUFactor::ETA ? "Eta" : "Forest-Tomlin")
             << " update"
             << std::endl;

   switch(pricing)
   {
   case 5 :
      pricer = new SPxWeightPR;
      break;
   case 4 :
      pricer = new SPxSteepPR;
      break;
   case 3 :
      pricer = new SPxHybridPR;
      break;
   case 2 :
      pricer = new SPxDevexPR;
      break;
   case 1 :
      pricer = new SPxParMultPR;
      break;
   case 0 : 
      /*FALLTHROUGH*/
   default :
      pricer = new SPxDefaultPR;
      break;
   }
   work.setPricer(pricer);

   std::cout << pricer->name() << " pricing" << std::endl;
   assert(work.isConsistent());

   switch(ratiotest)
   {
   case 2 :
      ratiotester = new SPxFastRT;
      std::cout << "Fast";
      break;
   case 1 :
      ratiotester = new SPxHarrisRT;
      std::cout << "Harris";
      break;
   case 0 :
      /*FALLTHROUGH*/
   default:
      ratiotester = new SPxDefaultRT;
      std::cout << "Default";
      break;
   }
   std::cout << " ratiotest" << std::endl;
   work.setTester(ratiotester);
   assert(work.isConsistent());

   switch(simplifing)
   {
   case 5 :
      simplifier = new SPxScale;
      std::cout << "Scale";
      break;
   case 4 : 
      simplifier = new SPxRedundantSM;
      std::cout << "Redundant";
      break;
   case 3 : 
      simplifier = new SPxRem1SM;
      std::cout << "Remove 1";
      break;
   case 2 :
      simplifier = new SPxAggregateSM;
      std::cout << "Aggregate";
      break;
   case 1 :
      simplifier = new SPxGeneralSM;
      std::cout << "General";
      break;
   case 0  :
      /*FALLTHROUGH*/
   default :
      std::cout << "No";
   }
   std::cout << " simplifier" << std::endl;
   work.setSimplifier(simplifier);
   assert(work.isConsistent());

   switch(starting)
   {
   case 3 :
      starter = new SPxVectorST;
      std::cout << "Vector";
      break;
   case 2 :
      starter = new SPxSumST;
      std::cout << "Sum";
      break;
   case 1 :
      starter = new SPxWeightST;
      std::cout << "Weight";
      break;
   case 0 :
      /*FALLTHROUGH*/
   default :
      std::cout << "No";
      break;
   }
   std::cout << " starter" << std::endl;
   work.setStarter(starter);
   assert(work.isConsistent());

   Timer timer;
   std::cout << "loading LP file " << filename << std::endl;

   if (!work.readFile(filename, &rownames, &colnames))
   {
      std::cout << "error while reading file" << std::endl;
      exit(1);
   }
   assert(work.isConsistent());

   std::cout << "LP has " 
             << work.nRows() 
             << "\trows and\n       "
             << work.nCols() 
             << "\tcolumns" 
             << std::endl;

   timer.start();
   std::cout << "solving LP" 
             << std::endl;

   work.solve();

   timer.stop();
   std::cout << "solution time  is: " 
             << timer.userTime() 
             << std::endl
             << "iterations    : " 
             << work.basis().iteration() 
             << std::endl;
   
   SoPlex::Status stat = work.status();

   switch (stat)
   {
   case SoPlex::SOLVED:
      std::cout << "solution value is: "
                << std::setprecision(10)
                << work.value();
      break;
   case SoPlex::UNBOUNDED:
      std::cout << "LP is unbounded";
      break;
   case SoPlex::INFEASIBLE:
      std::cout << "LP is infeasible";
      break;
   default:
      std::cout << "An error occurred during the solution process";
      break;
   }
   std::cout << std::endl;

   delete simplifier;
   delete starter;
   delete pricer;
   delete ratiotester;

   return 0;
}

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
