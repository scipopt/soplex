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
#pragma ident "@(#) $Id: example.cpp,v 1.70 2005/07/14 17:38:35 bzforlow Exp $"

#include <assert.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "spxdefines.h"
#include "soplex.h"
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
#include "spxintervalsm.h"
#include "spxaggregatesm.h"
#include "spxredundantsm.h"
#include "spxgeneralsm.h"
#include "spxscaler.h"
#include "spxequilisc.h"
#include "spxgeometsc.h"
#include "spxsumst.h"
#include "spxweightst.h"
#include "spxvectorst.h"
#include "slufactor.h"

using namespace soplex;

//------------------------------------------------------------------------
//    class MySoPlex
//------------------------------------------------------------------------

/** Here comes a simple derived class from #SoPlex.
 */
class MySoPlex : public SoPlex
{
public:
   /// default constructor
   MySoPlex(
      SPxSolver::Type p_type = SPxSolver::LEAVE, 
      SPxSolver::Representation p_rep = SPxSolver::COLUMN)
      : SoPlex(p_type, p_rep)
   {}
   //------------------------------------------------------------------------
   /// virtual destructor
   virtual ~MySoPlex()
   {}
   //------------------------------------------------------------------------
   void displayQuality() const
   {
      Real maxviol;
      Real sumviol;

      std::cout << "IEXAMP05 Violations (max/sum)" << std::endl;
                
      m_solver.qualConstraintViolation(maxviol, sumviol);

      std::cout << "IEXAMP06 Constraints      :" 
                << std::setw(16) << maxviol << "  " 
                << std::setw(16) << sumviol << std::endl;

      qualConstraintViolation(maxviol, sumviol);

      std::cout << "IEXAMP07       (unscaled) :" 
                << std::setw(16) << maxviol << "  " 
                << std::setw(16) << sumviol << std::endl;

      m_solver.qualBoundViolation(maxviol, sumviol);

      std::cout << "IEXAMP08 Bounds           :" 
                << std::setw(16) << maxviol << "  " 
                << std::setw(16) << sumviol << std::endl;

      qualBoundViolation(maxviol, sumviol);

      std::cout << "IEXAMP09       (unscaled) :" 
                << std::setw(16) << maxviol << "  " 
                << std::setw(16) << sumviol << std::endl;

      if (!m_vanished)
      {
         m_solver.qualSlackViolation(maxviol, sumviol);

         std::cout << "IEXAMP10 Slacks           :" 
                   << std::setw(16) << maxviol << "  " 
                   << std::setw(16) << sumviol << std::endl;

         m_solver.qualRedCostViolation(maxviol, sumviol);

         std::cout << "IEXAMP11 Reduced costs    :" 
                   << std::setw(16) << maxviol << "  " 
                   << std::setw(16) << sumviol << std::endl;
#if 0
         std::cout << "IEXAMP12 Proven dual bound:" 
                   << std::setw(20)
                   << std::setprecision(20)
                   << m_solver.provedDualbound() << std::endl;
#endif
      }
   }
   //------------------------------------------------------------------------
   void displayInfeasibility() const
   {
      assert(m_solver.status() == SPxSolver::INFEASIBLE);

#if 0
      if( m_solver.isProvenInfeasible() )
         std::cout << "IEXAMP13 Infeasibility is proven." << std::endl;
      else
         std::cout << "IEXAMP13 Infeasibility could not be proven!" << std::endl;
#endif
   }
};

//------------------------------------------------------------------------
//    main program
//------------------------------------------------------------------------

int main(int argc, const char* const argv[])
{
   const char* banner =
   "************************************************************************\n"
   "*                                                                      *\n"
   "*       SoPlex --- the Sequential object-oriented simPlex.             *\n"
   "*                  Release 1.2.2                                       *\n"
   "*    Copyright (C) 1997-1999 Roland Wunderling                         *\n"
   "*                  1997-2005 Konrad-Zuse-Zentrum                       *\n"
   "*                            fuer Informationstechnik Berlin           *\n"
   "*                                                                      *\n"
   "*  SoPlex is distributed under the terms of the ZIB Academic Licence.  *\n"
   "*  You should have received a copy of the ZIB Academic License         *\n"
   "*  along with SoPlex; If not email to soplex@zib.de.                   *\n"
   "*                                                                      *\n"
   "************************************************************************\n"
   ;

   const char* usage =
   "[options] LPfile [Basfile]\n\n"
   "          LPfile can be either in MPS or LPF format\n\n"
   "options:  (*) indicates default\n" 
   " -e        select entering algorithm (default is leaving)\n"
   " -r        select row wise representation (default is column)\n"
   " -i        select Eta-update (default is Forest-Tomlin)\n"
   " -x        output solution vector (works only together with -s0)\n"
   " -q        display solution quality\n"
   " -br       read file with starting basis from Basfile\n"
   " -bw       write file with optimal basis to Basfile\n"
   " -lSec     set timelimit to Sec seconds\n"
   " -dDelta   set maximal allowed bound violation to Delta\n"
   " -zzEps    set general zero tolerance to Eps\n\n"
   " -zfEps    set factorization zero tolerance to Eps\n\n"
   " -zuEps    set update zero tolerance to Eps\n\n"
   " -vLevel   set verbosity Level [0-3], default 1\n"
   " -V        show program version\n"
   " -h        show this help\n"
   "Simplifier:     Scaler:         Starter:     Pricer:        Ratiotester:\n"
   " -s0 none       -g0 none         -c0 none*   -p0 Textbook  -t0 Textbook\n"
   " -s1 General*   -g1 C-uni-Equi   -c1 Weight  -p1 ParMult   -t1 Harris\n"
   " -s2 Aggregate  -g2 R-uni-Equi   -c2 Sum     -p2 Devex     -t2 Fast*\n"
   " -s3 Redundant  -g3 bi-Equi      -c3 Vector  -p3 Hybrid\n"
   " -s4 Interval   -g4 bi-Equi+Geom1*           -p4 Steep*\n"
   "                -g5 bi-Equi+Geom8            -p5 Weight\n" 
   ;

   const char*               filename;
   char*                     basisname      = 0;
   SPxSolver::Type           type           = SPxSolver::LEAVE;
   SPxSolver::Representation representation = SPxSolver::COLUMN;
   SLUFactor::UpdateType     update         = SLUFactor::FOREST_TOMLIN;
   SPxSimplifier*            simplifier     = 0;
   SPxStarter*               starter        = 0;
   SPxPricer*                pricer         = 0;
   SPxRatioTester*           ratiotester    = 0;
   SPxScaler*                prescaler      = 0;
   SPxScaler*                postscaler     = 0;
   NameSet                   rownames;
   NameSet                   colnames;
   int                       starting       = 0;
   int                       pricing        = 4;
   int                       ratiotest      = 2;
   int                       scaling        = 4;
   int                       simplifing     = 1;
   Real                      timelimit      = -1.0;
   Real                      delta          = DEFAULT_BND_VIOL;
   Real                      epsilon        = DEFAULT_EPS_ZERO;
   Real                      epsilon_factor = DEFAULT_EPS_FACTOR;
   Real                      epsilon_update = DEFAULT_EPS_UPDATE;
   int                       verbose        = 1;
   bool                      print_solution = false;
   bool                      print_quality  = false;
   bool                      read_basis     = false;
   bool                      write_basis    = false;
   int                       precision;
   int                       optidx;

   for(optidx = 1; optidx < argc; optidx++)
   {
      if (*argv[optidx] != '-')
         break;

      switch(argv[optidx][1])
      {
      case 'b' :
         if (argv[optidx][2] == 'r')
            read_basis = true;
         if (argv[optidx][2] == 'w')
            write_basis = true;
         break;
      case 'c' :
         starting = atoi(&argv[optidx][2]);
         break;
      case 'd' :
         delta = atof(&argv[optidx][2]);
         break;
      case 'e':
         type = SPxSolver::ENTER;
         break;
      case 'g' :
         scaling = atoi(&argv[optidx][2]);
         break;
      case 'i' :
         update = SLUFactor::ETA;
         break;
      case 'l' :
         timelimit = atof(&argv[optidx][2]);
         break;
      case 'p' :
         pricing = atoi(&argv[optidx][2]);
         break;
      case 'q' :
         print_quality = true;
         break;
      case 'r' :
         representation = SPxSolver::ROW;
         break;
      case 's' :
         simplifing = atoi(&argv[optidx][2]);
         break;
      case 't' :
         ratiotest = atoi(&argv[optidx][2]);
         break;
      case 'v' :
         if (argv[optidx][2] >= '0' && argv[optidx][2] <= '9')
            verbose = argv[optidx][2] - '0';
         break;
      case 'V' :
         std::cout << banner << std::endl;
         exit(0);
      case 'x' :
         print_solution = true;
         break;
      case 'z' :
         switch(argv[optidx][2])
         {
         case 'z' :
            epsilon = atof(&argv[optidx][3]);
            break;
         case 'f' :
            epsilon_factor = atof(&argv[optidx][3]);
            break;
         case 'u' :
            epsilon_update = atof(&argv[optidx][3]);
            break;
         default :
            std::cout << "usage: " << argv[0] << " " << usage << std::endl;
            exit(0);
         }
         break;
      case 'h' :
      case '?' :
         std::cout << banner << std::endl;
         //lint -fallthrough
      default :
         std::cout << "usage: " << argv[0] << " " << usage << std::endl;
         exit(0);
      }
   }
   if ((argc - optidx) < 1 + (read_basis ? 1 : 0) + (write_basis ? 1 : 0))
   {
      std::cout << "usage: " << argv[0] << " " << usage << std::endl;
      exit(0);
   }
   filename  = argv[optidx];

   optidx++;

   if (read_basis || write_basis)
      basisname = strcpy(
         new char[strlen(argv[optidx]) + 1], argv[optidx]); 

   precision = int(-log10(delta)) + 1;

   Param::setEpsilon(epsilon);
   Param::setEpsilonFactorization(epsilon_factor);
   Param::setEpsilonUpdate(epsilon_update);
   Param::setVerbose(verbose);

   std::cout.setf(std::ios::scientific | std::ios::showpoint);
   std::cerr.setf(std::ios::scientific | std::ios::showpoint);

   // Send errors to a file (similarly for the other verbosity levels):
   // std::ofstream my_error_stream( "errors.txt" );
   // spxout.setStream( SPxOut::ERROR, my_error_stream );

   MySoPlex work(type, representation);

   work.setUtype(update);
   work.setTerminationTime(timelimit);
   work.setDelta(delta);

   std::cout << "IEXAMP12 Delta          = " << std::setw(16) << delta << std::endl
             << "IEXAMP13 Epsilon Zero   = " << std::setw(16) << Param::epsilon() << std::endl
             << "IEXAMP37 Epsilon Factor = " << std::setw(16) << Param::epsilonFactorization() << std::endl
             << "IEXAMP38 Epsilon Update = " << std::setw(16) << Param::epsilonUpdate() << std::endl;

   assert(work.isConsistent());

   std::cout << "IEXAMP14 "
             << (type == SPxSolver::ENTER ? "Entering" : "Leaving")
             << " algorithm" 
             << std::endl
             << "IEXAMP15 "
             << (representation == SPxSolver::ROW ? "Row" : "Column")
             << " representation" 
             << std::endl
             << "IEXAMP16 "
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

   std::cout << "IEXAMP17 " << pricer->getName() 
             << " pricing" << std::endl;

   assert(work.isConsistent());

   switch(ratiotest)
   {
   case 2 :
      ratiotester = new SPxFastRT;
      break;
   case 1 :
      ratiotester = new SPxHarrisRT;
      break;
   case 0 :
      /*FALLTHROUGH*/
   default:
      ratiotester = new SPxDefaultRT;
      break;
   }
   work.setTester(ratiotester);

   std::cout << "IEXAMP18 " << ratiotester->getName() 
             << " ratiotest" << std::endl;

   assert(work.isConsistent());

   switch(scaling)
   {
   case 5:
      prescaler  = new SPxEquiliSC(representation == SPxSolver::COLUMN, true);
      postscaler = new SPxGeometSC(representation == SPxSolver::COLUMN, 8);
      break; 
   case 4:
      prescaler  = new SPxEquiliSC(representation == SPxSolver::COLUMN, true);
      postscaler = new SPxGeometSC(representation == SPxSolver::COLUMN, 1);
      break; 
   case 3 :
      prescaler  = new SPxEquiliSC(representation == SPxSolver::COLUMN, true);
      postscaler = 0;
      break; 
   case 2 :
      prescaler  = new SPxEquiliSC(representation == SPxSolver::ROW, false);
      postscaler = 0;
      break; 
   case 1 :
      prescaler  = new SPxEquiliSC(representation == SPxSolver::COLUMN, false);
      postscaler = 0;
      break; 
   case 0 : 
      /*FALLTHROUGH*/
   default :
      prescaler  = 0;
      postscaler = 0;
      std::cout << "No";
      break;
   }
   work.setPreScaler(prescaler);
   work.setPostScaler(postscaler);

   std::cout << "IEXAMP19 "
             << ((prescaler != 0) ? prescaler->getName() : "no ") 
             << " / "
             << ((postscaler != 0) ? postscaler->getName() : "no ")
             << " scaling" << std::endl;

   assert(work.isConsistent());

   switch(simplifing)
   {
   case 4 :
      simplifier = new SPxIntervalSM;
      break;
   case 3 : 
      simplifier = new SPxRedundantSM;
      break;
   case 2 :
      simplifier = new SPxAggregateSM;
      break;
   case 1 :
      simplifier = new SPxGeneralSM;
      break;
   case 0  :
      /*FALLTHROUGH*/
   default :
      break;
   }
   work.setSimplifier(simplifier);

   std::cout << "IEXAMP20 "
             << ((simplifier == 0) ? "no" : simplifier->getName()) 
             << " simplifier" << std::endl;

   assert(work.isConsistent());

   switch(starting)
   {
   case 3 :
      starter = new SPxVectorST;
      break;
   case 2 :
      starter = new SPxSumST;
      break;
   case 1 :
      starter = new SPxWeightST;
      break;
   case 0 :
      /*FALLTHROUGH*/
   default :
      break;
   }
   work.setStarter(starter);

   std::cout << "IEXAMP21 "
             << ((starter == 0) ? "no" : starter->getName())
             << " starter" << std::endl;

   assert(work.isConsistent());

   Timer timer;
   std::cout << "IEXAMP22 loading LP file " << filename << std::endl;

   if (!work.readFile(filename, &rownames, &colnames))
   {
      ERROR( spxout << "EEXAMP23 error while reading file \"" 
                      << filename << "\"" << std::endl; )
      exit(1);
   }
   assert(work.isConsistent());

   std::cout << "IEXAMP24 LP has " 
             << work.nRows() << " rows "
             << work.nCols() << " columns " 
             << work.nNzos() << " nonzeros"
             << std::endl;

   // Should we read a basis ?
   if (read_basis)
   {
#if 0
      if (!work.readBasisFile(basisname, rownames, colnames))
      {
         ERROR( spxout << "EEXAMP25 error while reading file \"" 
                         << basisname << "\"" << std::endl; )
         exit(1);
      }
#endif
   }
   timer.start();
   std::cout << "IEXAMP26 solving LP" << std::endl;

   work.solve();

   timer.stop();

   std::cout << "IEXAMP01 Factorizations   : " << work.getFactorCount() << std::endl;
   std::cout << "IEXAMP02     Time spent   : " << work.getFactorTime() << std::endl;
   std::cout << "IEXAMP03 Solves           : " << work.getSolveCount() << std::endl;
   std::cout << "IEXAMP04     Time spent   : " << work.getSolveTime() << std::endl;

   std::cout << "IEXAMP27 solution time  is: " 
             << timer.userTime() 
             << std::endl
             << "IEXAMP28 iterations       : " 
             << work.iteration() 
             << std::endl;
   
   SPxSolver::Status stat = work.status();

   switch (stat)
   {
   case SPxSolver::OPTIMAL:
      std::cout << "IEXAMP29 solution value is: "
                << std::setprecision(precision)
                << work.objValue()
                << std::endl;

      if (print_quality)
         work.displayQuality();

      if (print_solution)
      {
         DVector objx(work.nCols());
            
         work.getPrimal(objx);
            
         for(int i = 0; i < work.nCols(); i++)
            if (isNotZero(objx[i], epsilon))
               std::cout << colnames[work.cId(i)] << "\t" 
                         << i << "\t"
                         << std::setw(16)
                         << std::setprecision(precision)
                         << objx[i] << std::endl;
         std::cout << "All other variable are zero." << std::endl;
      }
      if (write_basis)
         if (!work.writeBasisFile(basisname, rownames, colnames))
            ERROR( spxout << "EEXAMP30 error while writing file \"" 
                            << basisname << "\"" << std::endl; )
      break;
   case SPxSolver::UNBOUNDED:
      std::cout << "IEXAMP31 LP is unbounded" << std::endl;
      break;
   case SPxSolver::INFEASIBLE:
      std::cout << "IEXAMP32 LP is infeasible" << std::endl;
      if (print_quality)
         work.displayInfeasibility();
      break;
   case SPxSolver::ABORT_CYCLING:
      std::cout << "EEXAMP40 aborted due to cycling" << std::endl;
      break;
   case SPxSolver::ABORT_TIME:
      std::cout << "IEXAMP33 aborted due to time limit" << std::endl;
      break;
   case SPxSolver::ABORT_ITER:
      std::cout << "IEXAMP34 aborted due to iteration limit" << std::endl;
      break;
   case SPxSolver::ABORT_VALUE:
      std::cout << "IEXAMP35 aborted due to objective value limit" << std::endl;
      break;
   case SPxSolver::SINGULAR:
      ERROR( spxout << "EEXAMP39 basis is singular" << std::endl; )
      break;
   default:
      ERROR( spxout << "EEXAMP36 An error occurred during "
                      << "the solution process" << std::endl; )
      break;
   }
   std::cout << std::endl;

   if (prescaler != 0)
      delete prescaler;
   if (postscaler != 0)
      delete postscaler;
   if (simplifier != 0)
      delete simplifier;
   if (starter != 0)
      delete starter;

   assert(pricer != 0);
   delete pricer;

   assert(ratiotester != 0);
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
