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
#pragma ident "@(#) $Id: example.cpp,v 1.14 2002/01/11 21:05:30 bzfkocht Exp $"

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
/**@todo Flag to set epsilon and delta */
int main(int argc, char **argv)
{
   const char* usage = "[-eri][-bn][-ltime][-pn][-sn][-tn] mpsfile";
   const char*            filename;
   SoPlex::Type           type           = SoPlex::LEAVE;
   SoPlex::Representation representation = SoPlex::COLUMN;
   SLUFactor::UpdateType  update         = SLUFactor::FOREST_TOMLIN;
   int                    starter        = 0;
   int                    pricing        = 4;
   int                    ratiotest      = 2;
   int                    simplifier     = 3;
   double                 timelimit      = -1.0;
   int                    optind;

   for(optind = 1; optind < argc; optind++)
   {
      if (*argv[optind] != '-')
         break;

      switch(argv[optind][1])
      {
      case 'c' :
         starter = atoi(&argv[optind][2]);
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
         simplifier = atoi(&argv[optind][2]);
         break;
      case 't' :
         ratiotest = atoi(&argv[optind][2]);
         break;
      case '?' :
         std::cerr << "usage: " << argv[0] << ":" << usage << std::endl;
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

   switch (pricing)
   {
   case 5 :
      work.setPricer(new SPxWeightPR);
      std::cout << "Weight";
      break;
   case 4 :
      work.setPricer(new SPxSteepPR);
      std::cout << "Steepest edge";
      break;
   case 3 :
      work.setPricer(new SPxHybridPR);
      std::cout << "Hybrid";
      break;
   case 2 :
      work.setPricer(new SPxDevexPR);
      std::cout << "Devex";
      break;
   case 1 :
      work.setPricer(new SPxParMultPR);
      std::cout << "Partial multiple";
      break;
   case 0 :
      /*FALLTHROUGH*/
   default:
      work.setPricer(new SPxDefaultPR);
      std::cout << "Default";
      break;
   }
   std::cout << " pricing" << std::endl;
   assert(work.isConsistent());

   switch (ratiotest)
   {
   case 2 :
      work.setTester(new SPxFastRT);
      std::cout << "Fast";
      break;
   case 1 :
      work.setTester(new SPxHarrisRT);
      std::cout << "Harris";
      break;
   case 0 :
      /*FALLTHROUGH*/
   default:
      work.setTester(new SPxDefaultRT);
      std::cout << "Default";
      break;
   }
   std::cout << " ratiotest" << std::endl;
   assert(work.isConsistent());

   switch(simplifier)
   {
   case 5 :
      work.setSimplifier(new SPxScale);
      std::cout << "Scale";
      break;
   case 4 : 
      work.setSimplifier(new SPxRedundantSM);
      std::cout << "Redundant";
      break;
   case 3 : 
      work.setSimplifier(new SPxRem1SM);
      std::cout << "Remove 1";
      break;
   case 2 :
      work.setSimplifier(new SPxAggregateSM);
      std::cout << "Aggregate";
      break;
   case 1 :
      work.setSimplifier(new SPxGeneralSM);
      std::cout << "General";
      break;
   case 0  :
      /*FALLTHROUGH*/
   default :
      std::cout << "No";
   }
   std::cout << " simplifier" << std::endl;

   assert(work.isConsistent());

   switch(starter)
   {
   case 3 :
      work.setStarter(new SPxVectorST);
      std::cout << "Vector";
      break;
   case 2 :
      work.setStarter(new SPxSumST);
      std::cout << "Sum";
      break;
   case 1 :
      work.setStarter(new SPxWeightST);
      std::cout << "Weight";
      break;
   case 0 :
      /*FALLTHROUGH*/
   default :
      work.setStarter(static_cast<SPxStarter*>(0));
      std::cout << "Default";
      break;
   }

   std::cout << " starter" << std::endl;

   assert(work.isConsistent());

   Timer timer;
   std::cout << "loading LP file " << filename << std::endl;

   if (!work.readFile(filename))
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
