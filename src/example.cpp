/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: example.cpp,v 1.97 2007/09/04 10:24:10 bzfpfets Exp $"

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
#include "spxdantzigpr.h"
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
#include "spxmainsm.h"
#include "spxscaler.h"
#include "spxequilisc.h"
#include "spxgeometsc.h"
#include "spxsumst.h"
#include "spxweightst.h"
#include "spxvectorst.h"
#include "slufactor.h"
#include "spxout.h"

using namespace soplex;

//------------------------------------------------------------------------
//    class MySoPlex
//------------------------------------------------------------------------

/** LP solver example class for the command line.
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

      MSG_INFO1( spxout << "IEXAMP05 Violations (max/sum)" << std::endl; ) 
                
      m_solver.qualConstraintViolation(maxviol, sumviol);

      MSG_INFO1( spxout << "IEXAMP06 Constraints      :" 
                        << std::setw(16) << maxviol << "  " 
                        << std::setw(16) << sumviol << std::endl; )

      qualConstraintViolation(maxviol, sumviol);

      MSG_INFO1( spxout << "IEXAMP07       (unscaled) :" 
                        << std::setw(16) << maxviol << "  " 
                        << std::setw(16) << sumviol << std::endl; )

      m_solver.qualBoundViolation(maxviol, sumviol);

      MSG_INFO1( spxout << "IEXAMP08 Bounds           :" 
                        << std::setw(16) << maxviol << "  " 
                        << std::setw(16) << sumviol << std::endl; )

      qualBoundViolation(maxviol, sumviol);

      MSG_INFO1( spxout << "IEXAMP09       (unscaled) :" 
                        << std::setw(16) << maxviol << "  " 
                        << std::setw(16) << sumviol << std::endl; )

      if (!m_vanished)
      {
         m_solver.qualSlackViolation(maxviol, sumviol);

         MSG_INFO1( spxout << "IEXAMP10 Slacks           :" 
                           << std::setw(16) << maxviol << "  " 
                           << std::setw(16) << sumviol << std::endl; )

         m_solver.qualRedCostViolation(maxviol, sumviol);

         MSG_INFO1( spxout << "IEXAMP11 Reduced costs    :" 
                           << std::setw(16) << maxviol << "  " 
                           << std::setw(16) << sumviol << std::endl; )
#if 0
         MSG_INFO1( spxout << "IEXAMP12 Proven dual bound:" 
                           << std::setw(20)
                           << std::setprecision(20)
                           << m_solver.provedDualbound() << std::endl; )
#endif
      }
   }
   //------------------------------------------------------------------------
   void displayInfeasibility() const
   {
      assert(m_solver.status() == SPxSolver::INFEASIBLE);

#if 0
      if( m_solver.isProvenInfeasible() )
         MSG_INFO1( spxout << "IEXAMP13 Infeasibility is proven." 
                           << std::endl; )
      else
         MSG_INFO1( spxout << "IEXAMP13 Infeasibility could not be proven!"
                           << std::endl; )
#endif
   }
};


//------------------------------------------------------------------------
//    Helpers
//------------------------------------------------------------------------

void print_usage_and_exit( const char* const argv[] )
{
   /**
      The following code block is just to test compilation parameters.
   */
#if 0
   std::cout << "compiled with NDEBUG: "
   #ifdef NDEBUG
             << "yes"
   #else
             << "no"
   #endif
             << std::endl;

   std::cout << "compiled with WITH_WARNINGS: "
   #ifdef WITH_WARNINGS
             << "yes"
   #else
             << "no"
   #endif
             << std::endl;

   std::cout << "compiled with NO_ADDITIONAL_CHECKS: "
   #ifdef NO_ADDITIONAL_CHECKS
             << "yes (ENABLE_ADDITIONAL_CHECKS = " << ENABLE_ADDITIONAL_CHECKS << ")" 
   #else
             << "no (ENABLE_ADDITIONAL_CHECKS = " << ENABLE_ADDITIONAL_CHECKS << ")"
   #endif
             << std::endl;

   std::cout << "compiled with NO_CONSISTENCY_CHECKS: "
   #ifdef NO_CONSISTENCY_CHECKS
             << "yes"
   #else
             << "no"
   #endif
             << std::endl;

   std::cout << std::endl;
#endif

   const char* usage =
   "[options] LPfile [Basfile]\n\n"
   "          LPfile can be either in MPS or LPF format\n\n"
   "options:  (*) indicates default\n" 
   "          (!) indicates experimental features which may give wrong results\n" 
   " -e        select entering algorithm (default is leaving)\n"
   " -r        select row wise representation (default is column)\n"
   " -i        select Eta-update (default is Forest-Tomlin)\n"
   " -x        output solution vector (works only together with -s0)\n"
   " -q        display solution quality\n"
//   " -br       read file with starting basis from Basfile\n"
   " -bw       write file with optimal basis to Basfile\n"
   " -lSec     set timelimit to Sec seconds\n"
   " -dDelta   set maximal allowed bound violation to Delta\n"
   " -zzEps    set general zero tolerance to Eps\n\n"
   " -zfEps    set factorization zero tolerance to Eps\n\n"
   " -zuEps    set update zero tolerance to Eps\n\n"
   " -vLevel   set verbosity Level: from 0 (ERROR) to 5 (DEBUG), default 2\n"
   " -V        show program version\n"
   " -h        show this help\n"
   "Simplifier:     Scaler:         Starter:     Pricer:        Ratiotester:\n"
   " -s0 none       -g0 none         -c0 none*   -p0 Textbook  -t0 Textbook\n"
   " -s1 Main*      -g1 C-uni-Equi   -c1 Weight  -p1 ParMult   -t1 Harris\n"
   "                -g2 R-uni-Equi   -c2 Sum     -p2 Devex     -t2 Fast*\n"
   "                -g3 bi-Equi*     -c3 Vector  -p3 Hybrid!\n"
   "                -g4 bi-Equi+Geom1            -p4 Steep*\n"
   "                -g5 bi-Equi+Geom8            -p5 Weight\n"
   ;

   std::cerr << "usage: " << argv[0] << " " << usage << std::endl;
   exit(0);
}
//------------------------------------------------------------------------
void check_parameter( const char param, const char* const argv[] )
{
   if (param == '\0') 
      print_usage_and_exit( argv );
}
//------------------------------------------------------------------------
void print_algorithm_parameters(const MySoPlex&                 work,
                                const SPxSolver::Representation representation, 
                                const SLUFactor::UpdateType     update)
{
   MSG_INFO1( spxout 
                << "IEXAMP12 Delta          = " 
                << std::setw(16) << work.delta() << std::endl
                << "IEXAMP13 Epsilon Zero   = " 
                << std::setw(16) << Param::epsilon() << std::endl
                << "IEXAMP37 Epsilon Factor = " 
                << std::setw(16) << Param::epsilonFactorization() << std::endl
                << "IEXAMP38 Epsilon Update = " 
                << std::setw(16) << Param::epsilonUpdate() << std::endl
                << "IEXAMP14 "
                << (work.type() == SPxSolver::ENTER ? "Entering" : "Leaving")
                << " algorithm" << std::endl
                << "IEXAMP15 "
                << (representation == SPxSolver::ROW ? "Row" : "Column")
                << " representation" << std::endl
                << "IEXAMP16 "
                << (update == SLUFactor::ETA ? "Eta" : "Forest-Tomlin")
                << " update" << std::endl; )
}
//------------------------------------------------------------------------
void redirect_output( std::ostream&  myerrstream,
                      std::ostream&  myinfostream )
{
   myerrstream .setf( std::ios::scientific | std::ios::showpoint );
   myinfostream.setf( std::ios::scientific | std::ios::showpoint );
   spxout.setStream( SPxOut::ERROR,    myerrstream );
   spxout.setStream( SPxOut::WARNING,  myerrstream );
   spxout.setStream( SPxOut::INFO1,    myinfostream );
   spxout.setStream( SPxOut::INFO2,    myinfostream );
   spxout.setStream( SPxOut::INFO3,    myinfostream );
   spxout.setStream( SPxOut::DEBUG,    myinfostream );
}
//------------------------------------------------------------------------
void read_input_file(MySoPlex&      work,
                     const char*    filename,
                     NameSet&       rownames,
                     NameSet&       colnames)
{
   MSG_INFO1( spxout << "IEXAMP22 loading LP file " << filename << std::endl; )

   Timer timer;
   timer.start();

   if ( ! work.readFile(filename, &rownames, &colnames, NULL) )
   {
      MSG_INFO1( spxout << "EEXAMP23 error while reading file \"" 
                        << filename << "\"" << std::endl; )
      exit(1);
   }
   assert(work.isConsistent());

   timer.stop();

   MSG_INFO1( spxout << "IEXAMP24 LP has " 
                     << work.nRows() << " rows "
                     << work.nCols() << " columns " 
                     << work.nNzos() << " nonzeros"
                     << std::endl; )

   MSG_INFO1( spxout << "IEXAMP41 LP reading time: " << timer.userTime() << std::endl; )
}
//------------------------------------------------------------------------
void read_basis_file(MySoPlex&      /* work     */,
                     const char*    /* filename */,
                     NameSet&       /* rownames */,
                     NameSet&       /* colnames */)
{
#if 0
   if (!work.readBasisFile(basisname, rownames, colnames))
   {
      MSG_INFO1( spxout << "EEXAMP25 error while reading file \"" 
                        << basisname << "\"" << std::endl; )
      exit(1);
   }
#endif
}
//------------------------------------------------------------------------
SPxPricer* get_pricer(const int pricing)
{
   SPxPricer* pricer = NULL;
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
      pricer = new SPxDantzigPR;
      break;
   }

   assert(pricer != NULL);
   MSG_INFO1( spxout << "IEXAMP17 " << pricer->getName() 
                     << " pricing"  << std::endl; )
   return pricer;
}
//------------------------------------------------------------------------
SPxRatioTester* get_ratio_tester(const int ratiotest)
{
   SPxRatioTester* ratiotester = NULL;
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

   assert(ratiotester != NULL);
   MSG_INFO1( spxout << "IEXAMP18 " << ratiotester->getName() 
                     << " ratiotest" << std::endl; )
   return ratiotester;
}
//------------------------------------------------------------------------
void get_scalers(SPxScaler*& prescaler,
                 SPxScaler*& postscaler,
                 const int   scaling,
                 const int   representation)
{
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
      break;
   }

   MSG_INFO1( spxout << "IEXAMP19 "
                     << ((prescaler != 0) ? prescaler->getName() : "no ") 
                     << " / "
                     << ((postscaler != 0) ? postscaler->getName() : "no ")
                     << " scaling" << std::endl; )
}
//------------------------------------------------------------------------
SPxSimplifier* get_simplifier(const int simplifying)
{
   SPxSimplifier* simplifier = NULL;
   switch(simplifying)
   {
   case 1 :
      simplifier = new SPxMainSM;
      break;
   case 0  :
      /*FALLTHROUGH*/
   default :
      assert(simplifier == NULL);
      break;
   }

   MSG_INFO1( spxout << "IEXAMP20 "
                     << ((simplifier == 0) ? "no" : simplifier->getName()) 
                     << " simplifier" << std::endl; )
   return simplifier;
}
//------------------------------------------------------------------------
SPxStarter* get_starter(const int starting)
{
   SPxStarter* starter = NULL;
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

   MSG_INFO1( spxout << "IEXAMP21 "
                     << ((starter == 0) ? "no" : starter->getName())
                     << " starter" << std::endl; )

   return starter;
}
//------------------------------------------------------------------------
void solve_LP(MySoPlex&      work)
{
   Timer timer;
   timer.start();
   MSG_INFO1( spxout << "IEXAMP26 solving LP" << std::endl; )

   work.solve();
   timer.stop();

   MSG_INFO1( spxout 
      << "IEXAMP01 Factorizations   : " << work.getFactorCount() << std::endl
      << "IEXAMP02     Time spent   : " << work.getFactorTime() << std::endl
      << "IEXAMP03 Solves           : " << work.getSolveCount() << std::endl
      << "IEXAMP04     Time spent   : " << work.getSolveTime() << std::endl
      << "IEXAMP27 solution time  is: " << timer.userTime() << std::endl
      << "IEXAMP28 iterations       : " << work.iteration() << std::endl; )
}
//------------------------------------------------------------------------
void print_solution_and_status(const MySoPlex&      work,
                               const NameSet&       rownames,
                               const NameSet&       colnames,
                               const int            precision,
                               const bool           print_quality,
                               const bool           print_solution,
                               const bool           write_basis,
                               const char*          basisname)
{
   // get the solution status 
   SPxSolver::Status stat = work.status();

   switch (stat)
   {
   case SPxSolver::OPTIMAL:
      MSG_INFO1( spxout << "IEXAMP29 solution value is: "
                        << std::setprecision( precision )
                        << work.objValue() << std::endl; )

      if ( print_quality )
         work.displayQuality();

      if ( print_solution )
      {
         DVector objx(work.nCols());
            
         if( work.getPrimal(objx) != SPxSolver::ERROR )
         {
            for( int i = 0; i < work.nCols(); ++i ) 
            {
               if ( isNotZero( objx[i], Param::epsilon() ) )
                  MSG_INFO1( spxout << colnames[ work.cId(i) ] << "\t" 
                                    << i << "\t"
                                    << std::setw(16)
                                    << std::setprecision( precision )
                                    << objx[i] << std::endl; )
            }
            MSG_INFO1( spxout << "All other variable are zero." << std::endl; )
         }
      }
      if ( write_basis )
         if ( ! work.writeBasisFile( basisname, rownames, colnames ) )
            MSG_INFO1( spxout << "EEXAMP30 error while writing file \"" 
                              << basisname << "\"" << std::endl; )
      break;
   case SPxSolver::UNBOUNDED:
      MSG_INFO1( spxout << "IEXAMP31 LP is unbounded" << std::endl; )
      break;
   case SPxSolver::INFEASIBLE:
      MSG_INFO1( spxout << "IEXAMP32 LP is infeasible" << std::endl; )
      if ( print_solution )
      {
         DVector farkasx(work.nRows());
            
         if( work.getDualfarkas(farkasx) != SPxSolver::ERROR )
         {
            DVector proofvec(work.nCols());
            double lhs;
            double rhs;

            lhs = 0.0;
            rhs = 0.0;
            proofvec.clear();
            for( int i = 0; i < work.nRows(); ++i ) 
            {
               if ( isNotZero( farkasx[i], Param::epsilon() ) )
               {
                  MSG_INFO1( spxout << rownames[ work.rId(i) ] << "\t" 
                                    << i << "\t"
                                    << std::setw(16)
                                    << std::setprecision( precision )
                                    << farkasx[i] << "\t"; )
                  LPRow row;
                  work.getRow(i, row);
                  if( row.lhs() > -soplex::infinity )
                  {
                     MSG_INFO1( spxout << row.lhs() << " <= "; );
                  }
                  for( int j = 0; j < row.rowVector().size(); ++j )
                  {
                     if( row.rowVector().value(j) > 0 )
                        MSG_INFO1( spxout << "+"; );
                     MSG_INFO1( spxout 
                        << row.rowVector().value(j) << " "
                        << colnames[ work.cId(row.rowVector().index(j)) ] 
                        << " "; );
                  }
                  if( row.rhs() < soplex::infinity )
                  {
                     MSG_INFO1( spxout << "<= " << row.rhs(); );
                  }
                  MSG_INFO1( spxout << std::endl; );
                  if( farkasx[i] > 0.0 )
                  {
                     lhs += farkasx[i] * row.lhs();
                     rhs += farkasx[i] * row.rhs();
                  }
                  else
                  {
                     lhs += farkasx[i] * row.rhs();
                     rhs += farkasx[i] * row.lhs();
                  }
                  SVector vec(row.rowVector());
                  vec *= farkasx[i];
                  proofvec += vec;
               }
            }

            MSG_INFO1( spxout << "All other row multipliers are zero." << std::endl; )
            MSG_INFO1( spxout << "Farkas infeasibility proof: \t"; )
            MSG_INFO1( spxout << lhs << " <= "; )
            bool nonzerofound = false;
            for( int i = 0; i < work.nCols(); ++i )
            {
               if ( isNotZero( proofvec[i], Param::epsilon() ) )
               {
                  if( proofvec[i] > 0 )
                     MSG_INFO1( spxout << "+"; )
                  MSG_INFO1( spxout << proofvec[i] << " " 
                                        << colnames[ work.cId(i) ] << " "; )
                  nonzerofound = true;
               }
            }
            if( !nonzerofound )
               MSG_INFO1( spxout << "0 "; );
            MSG_INFO1( spxout << "<= " << rhs << std::endl; );
         }
      }
      if ( print_quality )
         work.displayInfeasibility();
      break;
   case SPxSolver::ABORT_CYCLING:
      MSG_INFO1( spxout << "EEXAMP40 aborted due to cycling" << std::endl; )
      break;
   case SPxSolver::ABORT_TIME:
      MSG_INFO1( spxout << "IEXAMP33 aborted due to time limit" << std::endl; )
      break;
   case SPxSolver::ABORT_ITER:
      MSG_INFO1( spxout << "IEXAMP34 aborted due to iteration limit" << std::endl; )
      break;
   case SPxSolver::ABORT_VALUE:
      MSG_INFO1( spxout << "IEXAMP35 aborted due to objective value limit" << std::endl; )
      break;
   case SPxSolver::SINGULAR:
      MSG_INFO1( spxout << "EEXAMP39 basis is singular" << std::endl; )
      break;
   default:
      MSG_INFO1( spxout << "EEXAMP36 An error occurred during "
                        << "the solution process" << std::endl; )
      break;
   }
   MSG_INFO1( spxout << std::endl; )
}
//------------------------------------------------------------------------
void clean_up( SPxScaler*&       prescaler, 
               SPxScaler*&       postscaler, 
               SPxSimplifier*&   simplifier, 
               SPxStarter*&      starter, 
               SPxPricer*&       pricer, 
               SPxRatioTester*&  ratiotester, 
               char*&            basisname )
{
   if ( prescaler != NULL ) {
      prescaler = NULL;
      delete prescaler;
   }
   if ( postscaler != NULL ) {
      postscaler = NULL;
      delete postscaler;
   }
   if ( simplifier != NULL ) {
      simplifier = NULL;
      delete simplifier;
   }
   if ( starter != 0 ) {
      starter = NULL;
      delete starter;
   }

   assert( pricer != 0 );
   pricer = NULL;
   delete pricer;

   assert( ratiotester != 0 );
   ratiotester = NULL;
   delete ratiotester;

   if ( basisname != 0 )
      delete [] basisname;
}


//------------------------------------------------------------------------
//    main program
//------------------------------------------------------------------------

int main(int argc, const char* const argv[])
{
   const char* banner1 =
   "************************************************************************\n"
   "*                                                                      *\n"
   "*       SoPlex --- the Sequential object-oriented simPlex.             *\n"
   "*                  Release ";
   const char* banner2 = "                                       *\n"
   "*    Copyright (C) 1997-1999 Roland Wunderling                         *\n"
   "*                  1997-2006 Konrad-Zuse-Zentrum                       *\n"
   "*                            fuer Informationstechnik Berlin           *\n"
   "*                                                                      *\n"
   "*  SoPlex is distributed under the terms of the ZIB Academic Licence.  *\n"
   "*  You should have received a copy of the ZIB Academic License         *\n"
   "*  along with SoPlex; If not email to soplex@zib.de.                   *\n"
   "*                                                                      *\n"
   "************************************************************************\n"
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
   int                       scaling        = 3;
   int                       simplifying    = 1;
   Real                      timelimit      = -1.0;
   Real                      delta          = DEFAULT_BND_VIOL;
   Real                      epsilon        = DEFAULT_EPS_ZERO;
   Real                      epsilon_factor = DEFAULT_EPS_FACTOR;
   Real                      epsilon_update = DEFAULT_EPS_UPDATE;
   int                       verbose        = SPxOut::INFO1;
   bool                      print_solution = false;
   bool                      print_quality  = false;
   bool                      read_basis     = false;
   bool                      write_basis    = false;
   int                       precision;
   int                       optidx;

   // parse the command line
   for(optidx = 1; optidx < argc; optidx++)
   {
      if (*argv[optidx] != '-')
         break;

      switch(argv[optidx][1])
      {
      case 'b' :
         check_parameter(argv[optidx][2], argv); // use -b{r,w}, not -b
//         if (argv[optidx][2] == 'r')
//            read_basis = true;   // not sure if reading a basis works
         if (argv[optidx][2] == 'w')
            write_basis = true;
         else
            print_usage_and_exit(argv);
         break;
      case 'c' :
         check_parameter(argv[optidx][2], argv); // use -c[0-3], not -c
         starting = atoi(&argv[optidx][2]);
         break;
      case 'd' :
         check_parameter(argv[optidx][2], argv); // use -dx, not -d
         delta = atof(&argv[optidx][2]);
         break;
      case 'e':
         type = SPxSolver::ENTER;
         break;
      case 'g' :
         check_parameter(argv[optidx][2], argv); // use -g[0-5], not -g
         scaling = atoi(&argv[optidx][2]);
         break;
      case 'i' :
         update = SLUFactor::ETA;
         break;
      case 'l' :
         if (argv[optidx][2] == '\0' )  // use -lx, not -l
            print_usage_and_exit( argv );
         timelimit = atof(&argv[optidx][2]);
         break;
      case 'p' :
         check_parameter(argv[optidx][2], argv); // use -p[0-5], not -p
         pricing = atoi(&argv[optidx][2]);
         break;
      case 'q' :
         print_quality = true;
         break;
      case 'r' :
         representation = SPxSolver::ROW;
         break;
      case 's' :
         check_parameter(argv[optidx][2], argv); // use -s[0-4], not -s
         simplifying = atoi(&argv[optidx][2]);
         break;
      case 't' :
         check_parameter(argv[optidx][2], argv); // use -r[0-2], not -r
         ratiotest = atoi(&argv[optidx][2]);
         break;
      case 'v' :
         check_parameter(argv[optidx][2], argv); // use -v[0-5], not -v
         if (argv[optidx][2] >= '0' && argv[optidx][2] <= '9')
            verbose = argv[optidx][2] - '0';
         break;
      case 'V' :
         std::cout << banner1 << SOPLEX_VERSION/100 << ".";
         std::cout << (SOPLEX_VERSION % 100)/10 << ".";
         std::cout << SOPLEX_VERSION % 10 << banner2 << std::endl;
         exit(0);
      case 'x' :
         print_solution = true;
         break;
      case 'z' :
         check_parameter(argv[optidx][2], argv); // must not be empty
         check_parameter(argv[optidx][3], argv); // must not be empty
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
            print_usage_and_exit( argv );
         }
         break;
      case 'h' :
      case '?' :
         std::cout << banner1 << SOPLEX_VERSION/100 << ".";
         std::cout << (SOPLEX_VERSION % 100)/10 << ".";
         std::cout << SOPLEX_VERSION % 10 << banner2 << std::endl;
         //lint -fallthrough
      default :
         print_usage_and_exit( argv );
      }
   }

   // enough arguments?
   if ((argc - optidx) < 1 + (read_basis ? 1 : 0) + (write_basis ? 1 : 0))
      print_usage_and_exit( argv );
   filename  = argv[optidx];

   ++optidx;

   if ( read_basis || write_basis )
      basisname = strcpy( new char[strlen(argv[optidx]) + 1], argv[optidx] ); 

   // Set some output parameters
   std::cout.setf( std::ios::scientific | std::ios::showpoint );
   std::cerr.setf( std::ios::scientific | std::ios::showpoint );
   precision = int(-log10(delta)) + 1;

#ifdef  SEND_ALL_OUTPUT_TO_FILES
   // Example of redirecting output to different files.
   // Default is cerr for errors and warnings, cout for everything else. 
   std::ofstream  myerrstream ( "errwarn.txt" );
   std::ofstream  myinfostream( "infos.txt" );
   redirect_output(myerrstream, myinfostream);
#endif

   // set some algorithm parameters
   Param::setEpsilon             ( epsilon );
   Param::setVerbose             ( verbose );
   Param::setEpsilonFactorization( epsilon_factor );
   Param::setEpsilonUpdate       ( epsilon_update );

   // create an instance of MySoPlex
   MySoPlex work(type, representation);
   work.setUtype                 ( update );
   work.setDelta                 ( delta  );
   work.setTerminationTime       ( timelimit );
   print_algorithm_parameters    ( work, representation, update );
   assert(work.isConsistent());

   // set pricer, starter, simplifier, and ratio tester
   work.setPricer    ( pricer      = get_pricer      (pricing) );
   work.setStarter   ( starter     = get_starter     (starting) );
   work.setSimplifier( simplifier  = get_simplifier  (simplifying) );
   work.setTester    ( ratiotester = get_ratio_tester(ratiotest) );
   assert(work.isConsistent());

   // set pre- and postscaler
   get_scalers(prescaler, postscaler, scaling, representation);
   work.setPreScaler (prescaler);
   work.setPostScaler(postscaler);
   assert(work.isConsistent());

   // read the LP from an input file (.lp or .mps)
   read_input_file(work, filename, rownames, colnames);

   // read a basis file if specified
   if (read_basis)
      read_basis_file(work, basisname, rownames, colnames);

   // solve the LP
   solve_LP(work);
   
   // print solution, status, infeasibility system,...
   print_solution_and_status(work, rownames, colnames, 
                             precision, print_quality,
                             print_solution, write_basis,
                             basisname);


   // clean up
   clean_up( prescaler, postscaler, simplifier, 
             starter, pricer, ratiotester, basisname );

   assert( pricer == 0 );

   return 0;
}

//---------------------------------------------------------------------
//  Example of constructing an LP using the callable library
//---------------------------------------------------------------------

#if 0
   SoPlex  mysoplex;
   NameSet rownames;
   NameSet colnames;

   const int num_rows = 5;
   const int num_cols = 5;

   // for all rows 
   for (int i = 0; i < num_rows; ++i) {

      // create a row and add the coefficient of each column
      DSVector dsvec(num_cols);
      for ( int j = 0; j < num_cols; ++j ) {
         const double mycoeff = (i+1) * (j+1);
         dsvec.add(j, mycoeff);
      }

      // create a row and add it to the LP
      const double myrhs = 10;
      LPRow row(-infinity, dsvec, myrhs);  // -infinity <= row <= myrhs
      mysoplex.addRow( row );

      // set the name of the row
      char myrowname[ 20 ];
      sprintf(myrowname, "row%d", i);
      rownames.add( myrowname );
   }

   // for all columns
   for (int j = 0; j < num_cols; ++j) {

      // set the name of the column
      char mycolname[ 20 ];
      sprintf(mycolname, "y%d", j);
      colnames.add( mycolname );

      // set the objective coefficient
      // (soplex always maximizes; multiply by -1 for minimizing)
      const double myobjcoeff = j+1;
      mysoplex.changeObj(j, myobjcoeff);
   }

   // bound the first variable (index 0) to 1 <= x <= 2
   mysoplex.changeLower(0, 1);
   mysoplex.changeUpper(0, 2);

   // dump the LP to stdout
   mysoplex.writeLPF(std::cout, &rownames, &colnames, NULL);
#endif

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
