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
#pragma ident "@(#) $Id: testcuts.cpp,v 1.1 2002/03/01 13:16:36 bzfpfend Exp $"

#include <assert.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "real.h"
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
      std::cout << iteration() << ":\t" << value() << std::endl;
      basis().desc().dump();
      return SoPlex::terminate();
   }

   void setUtype(SLUFactor::UpdateType tp)
   {
      m_slu.setUtype(tp);
   }
};

void printLine()
{
   std::cout << "------------------------------------------------------------"
             << std::endl;
}

int main(int argc, char **argv)
{
   const char* banner =
   "************************************************************************\n"
   "*                                                                      *\n"
   "*       SoPlex --- the Sequential object-oriented simPlex.             *\n"
   "*                  Release 1.2.1                                       *\n"
   "*    Copyright (C) 1997-1999 Roland Wunderling                         *\n"
   "*                  1997-2002 Konrad-Zuse-Zentrum                       *\n"
   "*                            fuer Informationstechnik Berlin           *\n"
   "*                                                                      *\n"
   "*  SoPlex is distributed under the terms of the ZIB Academic Licence.  *\n"
   "*  You should have received a copy of the ZIB Academic License         *\n"
   "*  along with SoPlex; If not email to soplex@zib.de.                   *\n"
   "*                                                                      *\n"
   "************************************************************************\n"
   ;

   const char* usage =
   "[options]\n\n"
   "options:  (*) indicates default\n" 
   " -e        select entering algorithm (default is leaving)\n"
   " -r        select row wise representation (default is column)\n"
   " -i        select Eta-update (default is Forest-Tomlin)\n"
   " -x        output solution vector (works only together with -s0)\n"
   " -bBasfile read file with starting basis\n"
   " -lSec     set timelimit to Sec seconds\n"
   " -dDelta   set maximal allowed bound violation to Delta\n"
   " -zZero    set zero tolerance to Zero\n\n"
   "Simplifier:         Starter:        Pricer:           Ratiotester:\n"
   " -s0  none           -c0  none*      -p0  Textbook     -t0  Textbook\n" 
   " -s1  General        -c1  Weight     -p1  ParMult      -t1  Harris\n"
   " -s2  Aggregate      -c2  Sum        -p2  Devex        -t2  Fast*\n"
   " -s3  Remove-1*      -c3  vector     -p3  Hybrid\n"
   " -s4  Redundant                      -p4  Steep*\n"
   " -s5  Scale                          -p5  Weight\n" 
   ;

   char*                  basisname      = 0;
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
   Real                   timelimit      = -1.0;
   Real                   delta          = DEFAULT_BND_VIOL;
   Real                   epsilon        = DEFAULT_EPS_ZERO;
   bool                   print_solution = false;
   int                    precision;
   int                    optind;

   for(optind = 1; optind < argc; optind++)
   {
      if (*argv[optind] != '-')
         break;

      switch(argv[optind][1])
      {
      case 'b' :
         basisname = strcpy(
            new char[strlen(&argv[optind][1]) + 1], &argv[optind][1]); 
         break;
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
      case 'x' :
         print_solution = true;
         break;
      case 'z' :
         epsilon = atof(&argv[optind][2]);
         break;
      case 'h' :
      case '?' :
         std::cout << banner << std::endl;
         /*FALLTHROUGH*/
      default :
         std::cerr << "usage: " << argv[0] << " " << usage << std::endl;
         exit(0);
      }
   }
   precision = int(-log10(delta)) + 1;

   Param::setEpsilon(epsilon);

   std::cout.setf(std::ios::scientific | std::ios::showpoint);

   MySoPlex work(type, representation);

   work.setUtype(update);
   work.setTerminationTime(timelimit);
   work.setDelta(delta);

   std::cout << "Delta   = " << std::setw(16) << delta << std::endl
             << "Epsilon = " << std::setw(16) 
             << Param::epsilon() << std::endl;

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

   /* Constructing test LP:
      min x1
      3 <=  x0 +  x1
      1 <= -x0 + 2x1
            x0 +  x1 <= 10
            x0 -  x1 <=  8
      0 <=  x0
      0 <=  x1

      valid cuts:  (a) 5 <= x0 + 2x1
                   (b) 2 <=       x1
   */
   std::cout << "initializing test LP" << std::endl;

   DSVector vec;
   
   // 3 <=  x0 +  x1
   vec.clear();
   vec.add( 0, +1.0 );
   vec.add( 1, +1.0 );
   work.addRow( LPRow( 3.0, vec, infinity ) );
   // 1 <= -x0 + 2x1
   vec.clear();
   vec.add( 0, -1.0 );
   vec.add( 1, +2.0 );
   work.addRow( LPRow( 1.0, vec, infinity ) );
   //       x0 +  x1 <= 10
   vec.clear();
   vec.add( 0, +1.0 );
   vec.add( 1, +1.0 );
   work.addRow( LPRow( -infinity, vec, 10.0 ) );
   //       x0 -  x1 <=  8
   vec.clear();
   vec.add( 0, +1.0 );
   vec.add( 1, -1.0 );
   work.addRow( LPRow( -infinity, vec, 8.0 ) );
   // bounds
   work.changeBounds( 0, 0.0, infinity );
   work.changeBounds( 1, 0.0, infinity );
   // objective
   work.changeObj( 0, 0.0 );
   work.changeObj( 1, 1.0 );
   work.changeSense( SoPlex::MINIMIZE );

   assert(work.isConsistent());

   std::cout << "LP has " 
             << work.nRows() 
             << "\trows and\n       "
             << work.nCols() 
             << "\tcolumns" 
             << std::endl;

   // Should we read a basis ?
   if (basisname != 0)
   {
      if (!work.readBasisFile(basisname, rownames, colnames))
      {
         std::cout << "error while reading file \"" 
                   << basisname << "\"" << std::endl;
         exit(1);
      }
   }


   SoPlex::Status stat;
   DVector        sol( 2 );
   Real           eps = 1e-4;
   std::cout << std::setprecision(2);

   /* solve LP */
   printLine();
   std::cout << "Solving..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::OPTIMAL );
   assert( EQ( work.value(), 4.0/3.0, eps ) );
   assert( EQ( sol[0], 5.0/3.0, eps ) );
   assert( EQ( sol[1], 4.0/3.0, eps ) );
   std::cout << std::endl;

   /* branching on x0:  x0 <= 1 */
   printLine();
   std::cout << "Branching on x0: x0 <= 1" << std::endl;
   work.changeUpper( 0, 1.0 );
   std::cout << "Resolving..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::OPTIMAL );
   assert( EQ( work.value(), 2.0, eps ) );
   assert( EQ( sol[0], 1.0, eps ) );
   assert( EQ( sol[1], 2.0, eps ) );
   std::cout << "Resetting upper bound on x0: 0 <= x0 <= infinity" << std::endl;
   work.changeUpper( 0, infinity );
   std::cout << std::endl;

   /* branching on x0:  x0 >= 2 */
   printLine();
   std::cout << "Branching on x0: x0 >= 2" << std::endl;
   work.changeLower( 0, 2.0 );
   std::cout << "Resolving..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::OPTIMAL );
   assert( EQ( work.value(), 1.5, eps ) );
   assert( EQ( sol[0], 2.0, eps ) );
   assert( EQ( sol[1], 1.5, eps ) );
   std::cout << "Resetting lower bound on x0: 0 <= x0 <= infinity" << std::endl;
   work.changeLower( 0, 0.0 );
   std::cout << std::endl;

   /* solve original LP again */
   printLine();
   std::cout << "Solving again original LP..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::OPTIMAL );
   assert( EQ( work.value(), 4.0/3.0, eps ) );
   assert( EQ( sol[0], 5.0/3.0, eps ) );
   assert( EQ( sol[1], 4.0/3.0, eps ) );
   std::cout << std::endl;

   /* fixing on x0:  x0 == 1 */
   printLine();
   std::cout << "Fixing on x0: x0 == 1" << std::endl;
   work.changeLower( 0, 1.0 );
   work.changeUpper( 0, 1.0 );
   std::cout << "Resolving..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::OPTIMAL );
   assert( EQ( work.value(), 2.0, eps ) );
   assert( EQ( sol[0], 1.0, eps ) );
   assert( EQ( sol[1], 2.0, eps ) );
   std::cout << "Resetting bounds on x0: 0 <= x0 <= infinity" << std::endl;
   work.changeLower( 0, 0.0 );
   work.changeUpper( 0, infinity );
   std::cout << std::endl;

   /* fixing on x0:  x0 >= 2 */
   printLine();
   std::cout << "Fixing on x0: x0 == 2" << std::endl;
   work.changeLower( 0, 2.0 );
   work.changeUpper( 0, 2.0 );
   std::cout << "Resolving..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::OPTIMAL );
   assert( EQ( work.value(), 1.5, eps ) );
   assert( EQ( sol[0], 2.0, eps ) );
   assert( EQ( sol[1], 1.5, eps ) );
   std::cout << "Resetting bounds on x0: 0 <= x0 <= infinity" << std::endl;
   work.changeLower( 0, 0.0 );
   work.changeUpper( 0, infinity );
   std::cout << std::endl;

   /* solve original LP again */
   printLine();
   std::cout << "Solving again original LP..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::OPTIMAL );
   assert( EQ( work.value(), 4.0/3.0, eps ) );
   assert( EQ( sol[0], 5.0/3.0, eps ) );
   assert( EQ( sol[1], 4.0/3.0, eps ) );
   std::cout << std::endl;

   /* branching on x1:  x1 <= 1 */
   printLine();
   std::cout << "Branching on x1: x0 <= 1" << std::endl;
   work.changeUpper( 1, 1.0 );
   std::cout << "Resolving..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::INFEASIBLE );
   std::cout << "Resetting upper bound on x1: 0 <= x0 <= infinity" << std::endl;
   work.changeUpper( 1, infinity );
   std::cout << std::endl;

   /* branching on x1:  x1 >= 2 */
   printLine();
   std::cout << "Branching on x1: x1 >= 2" << std::endl;
   work.changeLower( 1, 2.0 );
   std::cout << "Resolving..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::OPTIMAL );
   assert( EQ( work.value(), 2.0, eps ) );
   assert( EQ( sol[0], 1.0, eps ) ||
           EQ( sol[0], 3.0, eps ) );
   assert( EQ( sol[1], 2.0, eps ) );
   std::cout << "Resetting lower bound on x1: 0 <= x0 <= infinity" << std::endl;
   work.changeLower( 1, 0.0 );
   std::cout << std::endl;

   /* solve original LP again */
   printLine();
   std::cout << "Solving again original LP..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::OPTIMAL );
   assert( EQ( work.value(), 4.0/3.0, eps ) );
   assert( EQ( sol[0], 5.0/3.0, eps ) );
   assert( EQ( sol[1], 4.0/3.0, eps ) );
   std::cout << std::endl;

   /* fixing on x1:  x1 == 1 */
   printLine();
   std::cout << "Fixing on x1: x1 == 1" << std::endl;
   work.changeLower( 1, 1.0 );
   work.changeUpper( 1, 1.0 );
   std::cout << "Resolving..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::INFEASIBLE );
   std::cout << "Resetting bounds on x1: 0 <= x1 <= infinity" << std::endl;
   work.changeLower( 1, 0.0 );
   work.changeUpper( 1, infinity );
   std::cout << std::endl;

   /* fixing on x1:  x1 >= 2 */
   printLine();
   std::cout << "Fixing on x1: x1 == 2" << std::endl;
   work.changeLower( 1, 2.0 );
   work.changeUpper( 1, 2.0 );
   std::cout << "Resolving..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::OPTIMAL );
   assert( EQ( work.value(), 2.0, eps ) );
   assert( EQ( sol[0], 1.0, eps ) ||
           EQ( sol[0], 3.0, eps ) );
   assert( EQ( sol[1], 2.0, eps ) );
   std::cout << "Resetting bounds on x1: 0 <= x1 <= infinity" << std::endl;
   work.changeLower( 1, 0.0 );
   work.changeUpper( 1, infinity );
   std::cout << std::endl;

   /* solve original LP again */
   printLine();
   std::cout << "Solving again original LP..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::OPTIMAL );
   assert( EQ( work.value(), 4.0/3.0, eps ) );
   assert( EQ( sol[0], 5.0/3.0, eps ) );
   assert( EQ( sol[1], 4.0/3.0, eps ) );
   std::cout << std::endl;

   /* add cut 5 <= x0 + 2x1 */
   printLine();
   vec.clear();
   vec.add( 0, 1.0 );
   vec.add( 1, 2.0 );
   std::cout << "Adding cut 5 <= x0 + 2x1" << std::endl;
   work.addRow( LPRow( 5.0, vec, infinity ) );
   std::cout << "Resolving..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::OPTIMAL );
   assert( EQ( work.value(), 1.5, eps ) );
   assert( EQ( sol[0], 2.0, eps ) );
   assert( EQ( sol[1], 1.5, eps ) );
   std::cout << std::endl;

   /* add cut 2 <= x1 */
   printLine();
   vec.clear();
   vec.add( 1, 1.0 );
   std::cout << "Adding cut 2 <= x1" << std::endl;
   work.addRow( LPRow( 2.0, vec, infinity ) );
   std::cout << "Resolving..." << std::endl;
   std::cout << work << std::endl;
   work.solve();
   stat = work.status();
   std::cout << "Solution Status: " << int( stat ) << std::endl;
   work.basis().desc().dump();
   std::cout << "Objective value: " << work.value() << std::endl;
   work.getPrimal( sol );
   std::cout << "Solution       : " << sol[0] << ", " << sol[1] << std::endl;
   assert( stat == SoPlex::OPTIMAL );
   assert( EQ( work.value(), 2.0, eps ) );
   assert( EQ( sol[0], 1.0, eps ) ||
           EQ( sol[0], 3.0, eps ) );
   assert( EQ( sol[1], 2.0, eps ) );
   std::cout << std::endl;

   delete simplifier;
   delete starter;
   delete pricer;
   delete ratiotester;

   if (basisname != 0)
      delete [] basisname;

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
