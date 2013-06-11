/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <assert.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "spxdefines.h"
#include "soplex2.h"
#include "spxsolver.h"

#include "timer.h"
#include "spxgithash.h"
#include "spxpricer.h"
#include "spxdantzigpr.h"
#include "spxparmultpr.h"
#include "spxdevexpr.h"
#include "spxhybridpr.h"
#include "spxsteeppr.h"
#include "spxsteepexpr.h"
#include "spxweightpr.h"
#include "spxratiotester.h"
#include "spxharrisrt.h"
#include "spxdefaultrt.h"
#include "spxfastrt.h"
#include "spxboundflippingrt.h"
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
//    Helpers
//------------------------------------------------------------------------

static
void printVersionInfo( bool checkMode)
{
   const char* banner1 =
   "************************************************************************\n"
   "*                                                                      *\n"
   "*       SoPlex --- the Sequential object-oriented simPlex.             *\n"
   ;

   const char* banner2 =
   "*                                                                      *\n"
   "*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                       *\n"
   "*                            fuer Informationstechnik Berlin           *\n"
   "*                                                                      *\n"
   "*  SoPlex is distributed under the terms of the ZIB Academic Licence.  *\n"
   "*  You should have received a copy of the ZIB Academic License         *\n"
   "*  along with SoPlex; If not email to soplex@zib.de.                   *\n"
   "*                                                                      *\n"
   "************************************************************************\n"
   ;

   if( !checkMode )
      std::cout << banner1;

#if (SOPLEX_SUBVERSION > 0)
   if( !checkMode )
      std::cout <<    "*                  Version ";
   else
      std::cout << "SoPlex version ";
   std::cout << SOPLEX_VERSION/100 << "."
             << (SOPLEX_VERSION % 100)/10 << "."
             << SOPLEX_VERSION % 10 << "."
             << SOPLEX_SUBVERSION
             << " - Githash "
             << std::setw(13) << std::setiosflags(std::ios::left) << getGitHash();
   if( !checkMode )
      std::cout << "             *\n" << banner2 << std::endl;
   else
      std::cout << "\n";
#else
   if( !checkMode )
      std::cout <<    "*                  Release ";
   else
      std::cout << "SoPlex release ";
   std::cout << SOPLEX_VERSION/100 << "."
             << (SOPLEX_VERSION % 100)/10 << "."
             << SOPLEX_VERSION % 10
             << " - Githash "
             << std::setw(13) << std::setiosflags(std::ios::left) << getGitHash();
   if( !checkMode )
      std::cout << "               *\n" << banner2 << std::endl;
   else
      std::cout << "\n";
#endif

   /// The following code block is tests and shows compilation parameters.
   std::cout << "[NDEBUG:"
#ifdef NDEBUG
             << "YES"
#else
             << "NO"
#endif
             << "]";

   std::cout << "[WITH_WARNINGS:"
#ifdef WITH_WARNINGS
             << "YES"
#else
             << "NO"
#endif
             << "]";

   std::cout << "[ENABLE_ADDITIONAL_CHECKS:"
#ifdef ENABLE_ADDITIONAL_CHECKS
             << "YES"
#else
             << "NO"
#endif
             << "]";

   std::cout << "[ENABLE_CONSISTENCY_CHECKS:"
#ifdef ENABLE_CONSISTENCY_CHECKS
             << "YES"
#else
             << "NO"
#endif
             << "]";

   std::cout << "[SOPLEX_WITH_GMP:"
#ifdef SOPLEX_WITH_GMP
             << "YES"
#else
             << "NO"
#endif
             << "]" << std::endl;

   std::cout << std::endl;
}

#if 0
static
void print_short_version_info()
{
   const char* banner1 =
   "************************************************************************\n"
   "* SoPlex --- the Sequential object-oriented simPlex. ";
   const char* banner2 =
   "* Copyright (C)  1996-2013 Zuse Institute Berlin                       *\n"
   "************************************************************************\n";

   std::cout << banner1;
#if (SOPLEX_SUBVERSION > 0)
   std::cout <<    "Version "
             << SOPLEX_VERSION/100 << "."
             << (SOPLEX_VERSION % 100)/10 << "."
             << SOPLEX_VERSION % 10 << "."
             << SOPLEX_SUBVERSION
             << "   *\n";
#else
   std::cout <<    "Release "
             << SOPLEX_VERSION/100 << "."
             << (SOPLEX_VERSION % 100)/10 << "."
             << SOPLEX_VERSION % 10
             << "     *\n";
#endif
   std::cout << banner2 << std::endl;
}
#endif

//------------------------------------------------------------------------
static
void printUsage( const char* const argv[] )
{
   const char* usage =
      "[options] LPfile [Basfile]\n\n"
      "          LPfile can be either in MPS or LPF format\n\n"
      "options:  (*) indicates default\n"
      "          (!) indicates experimental features which may give wrong results\n"
      " -X        solve rational LP\n"
      " -e        select entering algorithm (default is leaving)\n"
      " -r        select row wise representation (default is column)\n"
      " -i        select Eta-update (default is Forest-Tomlin)\n"
      " -x        output solution vector\n"
      " -y        output dual multipliers\n"
      " -q        display solution quality\n"
      " -br       read file with starting basis from Basfile\n"
      " -bw       write file with optimal basis to Basfile\n"
      " -l        set time limit in seconds\n"
      " -L        set iteration limit\n"
      " -f        set primal feasibility tolerance\n"
      " -o        set optimality, i.e., dual feasibility tolerance\n"
      " -d        set primal and dual feasibility tolerance to same value\n"
      " -R        set threshold for tolerances below which iterative refinement is applied\n"
      " -zz       set general zero tolerance\n"
      " -zf       set factorization zero tolerance\n"
      " -zu       set update zero tolerance\n"
      " -P        enable partial (= incomplete) pricing for leaving algorithm\n"
      " -v        set verbosity Level: from 0 (ERROR) to 5 (INFO3), default 3 (INFO1)\n"
      " -V        show program version\n"
      " -C        check mode (for check scripts)\n"
      " -h        show this help\n\n"
      "Simplifier:  Scaler:           Starter:    Pricer:        Ratiotester:\n"
      " -s0 none     -g0 none          -c0 none*   -p0 Textbook   -t0 Textbook\n"
      " -s1 Main*    -g1 uni-Equi      -c1 Weight  -p1 ParMult    -t1 Harris\n"
      "              -g2 bi-Equi*      -c2 Sum     -p2 Devex      -t2 Fast*\n"
      "              -g3 bi-Equi+Geo1  -c3 Vector  -p3 Hybrid!    -t3 Bound Flipping\n"
      "              -g4 bi-Equi+Geo8              -p4 Steep*\n"
      "                                            -p5 Weight\n"
      "                                            -p6 SteepExactSetup\n"
      ;

   std::cerr << "usage: " << argv[0] << " " << usage << std::endl;
   exit(0);
}

static
void displayQuality( SoPlex2& SoPlexShell, bool checkMode )
{

   /// @todo needs more functionality from soplex2 regarding original/simplified LP
   #if 0
   Real maxviol;
   Real sumviol;

   if ( checkMode )
   {
      MSG_INFO1( spxout << "IEXAMP05 Violations (max/sum)" << std::endl; )

      // get violations of simpified LP
      SoPlexShell.getConstraintViolationReal(primal, maxviol, sumviol);
      m_solver.qualConstraintViolation(maxviol, sumviol);

      MSG_INFO1( spxout << "IEXAMP06 Constraints      :"
      << std::setw(16) << maxviol << "  "
      << std::setw(16) << sumviol << std::endl; )

      // get violations of original LP
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
else
{
   MSG_INFO1( spxout << "Violations (max/sum)" << std::endl; )

   m_solver.qualConstraintViolation(maxviol, sumviol);

   MSG_INFO1( spxout << "Constraints      :"
   << std::setw(16) << maxviol << "  "
   << std::setw(16) << sumviol << std::endl; )

   qualConstraintViolation(maxviol, sumviol);

   MSG_INFO1( spxout << "      (unscaled) :"
   << std::setw(16) << maxviol << "  "
   << std::setw(16) << sumviol << std::endl; )

   m_solver.qualBoundViolation(maxviol, sumviol);

   MSG_INFO1( spxout << "Bounds           :"
   << std::setw(16) << maxviol << "  "
   << std::setw(16) << sumviol << std::endl; )

   qualBoundViolation(maxviol, sumviol);

   MSG_INFO1( spxout << "      (unscaled) :"
   << std::setw(16) << maxviol << "  "
   << std::setw(16) << sumviol << std::endl; )

   if (!m_vanished)
   {
      m_solver.qualSlackViolation(maxviol, sumviol);

      MSG_INFO1( spxout << "Slacks           :"
      << std::setw(16) << maxviol << "  "
      << std::setw(16) << sumviol << std::endl; )

      m_solver.qualRedCostViolation(maxviol, sumviol);

      MSG_INFO1( spxout << "Reduced costs    :"
      << std::setw(16) << maxviol << "  "
      << std::setw(16) << sumviol << std::endl; )
      #if 0
      MSG_INFO1( spxout << "Proven dual bound:"
      << std::setw(20)
      << std::setprecision(20)
      << m_solver.provedDualbound() << std::endl; )
      #endif
}
}
#endif
}

//------------------------------------------------------------------------
static
void displayInfeasibility(SoPlex2& SoPlexShell)
{
   assert(SoPlexShell.statusReal() == SPxSolver::INFEASIBLE);

   #if 0
   if ( checkMode )
   {
      if( m_solver.isProvenInfeasible() )
         MSG_INFO1( spxout << "IEXAMP13 Infeasibility is proven." << std::endl; )
         else
            MSG_INFO1( spxout << "IEXAMP13 Infeasibility could not be proven!" << std::endl; )
   }
   else
   {
      if ( m_solver.isProvenInfeasible() )
      {
         MSG_INFO1( spxout << "Infeasibility is proven." << std::endl; )
      }
      else
      {
         MSG_INFO1( spxout << "Infeasibility could not be proven!" << std::endl; )
      }
   }
   #endif
}


//------------------------------------------------------------------------
static
void checkParameter(const char param, const char* const argv[])
{
   if (param == '\0')
      printUsage( argv );
}

//------------------------------------------------------------------------
static
void printAlgorithmParameters(
   SoPlex2& SoPlexShell,
   bool checkMode)
{
   if( checkMode )
   {
      MSG_INFO1( spxout
    << "IEXAMP12 Feastol        = "
    << std::setw(16) << SoPlexShell.rationalParam(SoPlex2::FEASTOL) << std::endl
    << "IEXAMP52 Opttol         = "
    << std::setw(16) << SoPlexShell.rationalParam(SoPlex2::OPTTOL) << std::endl
    << "IEXAMP53 Irthreshold    = "
    << std::setw(16) << SoPlexShell.realParam(SoPlex2::IRTHRESHOLD) << std::endl
    << "IEXAMP13 Epsilon Zero   = "
    << std::setw(16) << SoPlexShell.realParam(SoPlex2::EPSILON_ZERO) << std::endl
    << "IEXAMP37 Epsilon Factor = "
    << std::setw(16) << SoPlexShell.realParam(SoPlex2::EPSILON_FACTORIZATION) << std::endl
    << "IEXAMP38 Epsilon Update = "
    << std::setw(16) << SoPlexShell.realParam(SoPlex2::EPSILON_UPDATE) << std::endl
    << "IEXAMP14 "
    << (SoPlexShell.intParam(SoPlex2::ALGORITHM) == SoPlex2::ALGORITHM_ENTER ? "Entering" : "Leaving")
    << " algorithm" << std::endl
    << "IEXAMP15 "
    << (SoPlexShell.intParam(SoPlex2::REPRESENTATION) == SoPlex2::REPRESENTATION_ROW ? "Row" : "Column")
    << " representation" << std::endl
    << "IEXAMP16 "
    << (SoPlexShell.intParam(SoPlex2::FACTOR_UPDATE_TYPE) == SoPlex2::FACTOR_UPDATE_TYPE_ETA ? "Eta" : "Forest-Tomlin")
    << " update" << std::endl; )
   }
   else
   {
      MSG_INFO1( spxout
    << "SoPlex parameters: " << std::endl
    << "Feastol        = "
    << std::setw(16) << SoPlexShell.rationalParam(SoPlex2::FEASTOL) << std::endl
    << "Opttol         = "
    << std::setw(16) << SoPlexShell.rationalParam(SoPlex2::OPTTOL) << std::endl
    << "Irthreshold    = "
    << std::setw(16) << SoPlexShell.realParam(SoPlex2::IRTHRESHOLD) << std::endl
    << "Epsilon Zero   = "
    << std::setw(16) << SoPlexShell.realParam(SoPlex2::EPSILON_ZERO) << std::endl
    << "Epsilon Factor = "
    << std::setw(16) << SoPlexShell.realParam(SoPlex2::EPSILON_FACTORIZATION) << std::endl
    << "Epsilon Update = "
    << std::setw(16) << SoPlexShell.realParam(SoPlex2::EPSILON_UPDATE) << std::endl
    << std::endl
    << "algorithm      = " << (SoPlexShell.intParam(SoPlex2::ALGORITHM) == SoPlex2::ALGORITHM_ENTER ? "Entering" : "Leaving")
    << std::endl
    << "representation = " << (SoPlexShell.intParam(SoPlex2::REPRESENTATION) == SoPlex2::REPRESENTATION_ROW ? "Row" : "Column")
    << std::endl
    << "update         = " << (SoPlexShell.intParam(SoPlex2::FACTOR_UPDATE_TYPE) == SoPlex2::FACTOR_UPDATE_TYPE_ETA ? "Eta" : "Forest-Tomlin")
    << std::endl; )
   }

   // @todo do we need checkMode keys for these lines?
   MSG_INFO1( spxout
    << "pricer         = " << SoPlexShell.getPricerName()
    << std::endl
    << "starter        = " << SoPlexShell.getStarterName()
    << std::endl
    << "simplifier     = " << SoPlexShell.getSimplifierName()
    << std::endl
    << "ratiotest      = " << SoPlexShell.getRatiotesterName()
    << std::endl
    << "scaling        = " << SoPlexShell.getFirstScalerName() << " / " << SoPlexShell.getSecondScalerName()
    << std::endl
    << std::endl; )
}

//------------------------------------------------------------------------
static
void setPricer(
   SoPlex2& SoPlexShell,
   const int pricing)
{
   /// @todo weight and partial pricers are not available
   switch(pricing)
   {
   case 6 :
      SoPlexShell.setIntParam(SoPlex2::PRICER, SoPlex2::PRICER_STEEP);
      break;
   case 5 :
      SoPlexShell.setIntParam(SoPlex2::PRICER, SoPlex2::PRICER_AUTO);
      break;
   case 4 :
      SoPlexShell.setIntParam(SoPlex2::PRICER, SoPlex2::PRICER_QUICKSTEEP);
      break;
   case 3 :
      SoPlexShell.setIntParam(SoPlex2::PRICER, SoPlex2::PRICER_HYBRID);
      break;
   case 2 :
      SoPlexShell.setIntParam(SoPlex2::PRICER, SoPlex2::PRICER_DEVEX);
      break;
   case 1 :
      SoPlexShell.setIntParam(SoPlex2::PRICER, SoPlex2::PRICER_PARMULT);
      break;
   case 0 :
      SoPlexShell.setIntParam(SoPlex2::PRICER, SoPlex2::PRICER_DANTZIG);
      break;
      /*FALLTHROUGH*/
   default :
      SoPlexShell.setIntParam(SoPlex2::PRICER, SoPlex2::PRICER_AUTO);
      break;
   }
}

//------------------------------------------------------------------------
static
void setRatiotester(
   SoPlex2& SoPlexShell,
   const int ratiotest)
{
   switch(ratiotest)
   {
   case 3 :
      SoPlexShell.setIntParam(SoPlex2::RATIOTESTER, SoPlex2::RATIOTESTER_BOUNDFLIPPING);
      break;
   case 2 :
      SoPlexShell.setIntParam(SoPlex2::RATIOTESTER, SoPlex2::RATIOTESTER_FAST);
      break;
   case 1 :
      SoPlexShell.setIntParam(SoPlex2::RATIOTESTER, SoPlex2::RATIOTESTER_HARRIS);
      break;
   case 0 :
      /*FALLTHROUGH*/
   default:
      SoPlexShell.setIntParam(SoPlex2::RATIOTESTER, SoPlex2::RATIOTESTER_TEXTBOOK);
      break;
   }
}

//------------------------------------------------------------------------
static
void setScalers(
   SoPlex2&    SoPlexShell,
   const int   scaling)
{
   switch(scaling)
   {
   case 4:
      SoPlexShell.setIntParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER, SoPlex2::SCALER_BIEQUI);
      SoPlexShell.setIntParam(SoPlex2::SCALER_AFTER_SIMPLIFIER, SoPlex2::SCALER_GEO8);
      break;
   case 3:
      SoPlexShell.setIntParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER, SoPlex2::SCALER_BIEQUI);
      SoPlexShell.setIntParam(SoPlex2::SCALER_AFTER_SIMPLIFIER, SoPlex2::SCALER_GEO1);
      break;
   case 2 :
      SoPlexShell.setIntParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER, SoPlex2::SCALER_BIEQUI);
      SoPlexShell.setIntParam(SoPlex2::SCALER_AFTER_SIMPLIFIER, SoPlex2::SCALER_OFF);
      break;
   case 1 :
      SoPlexShell.setIntParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER, SoPlex2::SCALER_UNIEQUI);
      SoPlexShell.setIntParam(SoPlex2::SCALER_AFTER_SIMPLIFIER, SoPlex2::SCALER_OFF);
      break;
   case 0 :
      /*FALLTHROUGH*/
   default :
      SoPlexShell.setIntParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER, SoPlex2::SCALER_OFF);
      SoPlexShell.setIntParam(SoPlex2::SCALER_AFTER_SIMPLIFIER, SoPlex2::SCALER_OFF);
      break;
   }
}

//------------------------------------------------------------------------
static
void setStarter(
   SoPlex2& SoPlexShell,
   const int starting)
{
   switch(starting)
   {
   case 3 :
      SoPlexShell.setIntParam(SoPlex2::STARTER, SoPlex2::STARTER_VECTOR);
      break;
   case 2 :
      SoPlexShell.setIntParam(SoPlex2::STARTER, SoPlex2::STARTER_SUM);
      break;
   case 1 :
      SoPlexShell.setIntParam(SoPlex2::STARTER, SoPlex2::STARTER_WEIGHT);
      break;
   case 0 :
      /*FALLTHROUGH*/
   default :
      SoPlexShell.setIntParam(SoPlex2::STARTER, SoPlex2::STARTER_OFF);
      break;
   }
}

//------------------------------------------------------------------------
#ifdef SEND_ALL_OUTPUT_TO_FILES
static
void redirectOutput(
   std::ostream&  myerrstream,
   std::ostream&  myinfostream
   )
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
#endif

//------------------------------------------------------------------------
static
void printSolutionAndStatus(
   SoPlex2&             SoPlexShell,
   const NameSet&       rownames,
   const NameSet&       colnames,
   const int            precision,
   const bool           print_quality,
   const bool           print_solution,
   const bool           print_dual,
   const bool           write_basis,
   const char*          basisname,
   bool                 checkMode)
{
   // get the solution status
   SPxSolver::Status status = SoPlexShell.statusReal();

   if( !checkMode )
      MSG_INFO1( spxout << std::endl; )
   switch (status)
   {
   case SPxSolver::OPTIMAL:
      if( checkMode )
    MSG_INFO1( spxout << "IEXAMP29 solution value is: " << std::setprecision( precision ) << SoPlexShell.objValueReal() << std::endl; )
      else
         MSG_INFO1( spxout << "Solution value is: " << std::setprecision( precision ) << SoPlexShell.objValueReal() << std::endl; )

      if( print_quality )
         displayQuality( SoPlexShell, checkMode );

      if( print_solution )
      {
         DVector objx(SoPlexShell.numColsReal());

         if( SoPlexShell.getPrimalReal(objx) )
         {
            MSG_INFO1( spxout << std::endl << "Primal solution (name, id, value):" << std::endl; )
            for( int i = 0; i < SoPlexShell.numColsReal(); ++i )
            {
               if( isNotZero( objx[i], 0.001 * Real(SoPlexShell.rationalParam(SoPlex2::FEASTOL)) ) )
                  MSG_INFO1( spxout << colnames[ SoPlexShell.colIdReal(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(17)
                                    << std::setprecision( precision )
                                    << objx[i] << std::endl; )
            }
            MSG_INFO1( spxout << "All other variables are zero (within " << std::setprecision(1) << 0.001 * SoPlexShell.rationalParam(SoPlex2::FEASTOL) << ")." << std::endl; )
         }
      }
      if( print_dual )
      {
         DVector objy( SoPlexShell.numRowsReal() );
         bool allzero = true;

         if( SoPlexShell.getDualReal(objy) )
         {
            MSG_INFO1( spxout << std::endl << "Dual multipliers (name, id, value):" << std::endl; )
            for( int i = 0; i < SoPlexShell.numRowsReal(); ++i )
            {
               if( isNotZero(objy[i] , 0.001 * Real(SoPlexShell.rationalParam(SoPlex2::OPTTOL))) )
               {
                  MSG_INFO1( spxout << rownames[ SoPlexShell.rowIdReal(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(17)
                                    << std::setprecision( precision )
                                    << objy[i] << std::endl; )
                  allzero = false;
               }
            }

            MSG_INFO1( spxout << "All " << (allzero ? "" : "other ") << "dual values are zero (within "
            << std::setprecision(1) << 0.001 * SoPlexShell.rationalParam(SoPlex2::OPTTOL) << ")." << std::endl; )

            if( !allzero )
            {
               if( SoPlexShell.intParam(SoPlex2::OBJSENSE) == SoPlex2::OBJSENSE_MINIMIZE )
               {
                  MSG_INFO1( spxout << "Minimizing: a positive/negative value corresponds to left-hand (>=) resp. right-hand (<=) side."
                                    << std::endl; )
               }
               else
               {
                  MSG_INFO1( spxout << "Maximizing: a positive/negative value corresponds to right-hand (<=) resp. left-hand (>=) side."
                                    << std::endl; )
               }
            }
         }
      }
      if( write_basis )
      {
         MSG_INFO1( spxout << "Writing basis of original problem to file " << basisname << std::endl; )
         SoPlexShell.writeBasisFileReal( basisname, &rownames, &colnames );
      }
      break;
   case SPxSolver::UNBOUNDED:
      if( checkMode )
    MSG_INFO1( spxout << "IEXAMP31 LP is unbounded" << std::endl; )
      else
    MSG_INFO1( spxout << "LP is unbounded" << std::endl; )

      if( print_solution )
      {
         DVector objx(SoPlexShell.numColsReal());
         if( SoPlexShell.getPrimalReal(objx) )
         {
            MSG_INFO1( spxout << std::endl << "Primal solution (name, id, value):" << std::endl; )
            for( int i = 0; i < SoPlexShell.numColsReal(); ++i )
            {
               if( isNotZero( objx[i], 0.001 * Real(SoPlexShell.rationalParam(SoPlex2::FEASTOL)) ) )
                  MSG_INFO1( spxout << colnames[ SoPlexShell.colIdReal(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(17)
                                    << std::setprecision( precision )
                                    << objx[i] << std::endl; )
            }
            MSG_INFO1( spxout << "All other variables are zero (within " << std::setprecision(1) << 0.001 * SoPlexShell.rationalParam(SoPlex2::FEASTOL) << ")." << std::endl; )
         }

         DVector objcoef(SoPlexShell.numColsReal());
         DVector ray(SoPlexShell.numColsReal());
         if( SoPlexShell.getPrimalrayReal(ray) )
         {
            Real rayobjval = 0.0;

            SoPlexShell.getObjReal(objcoef);

            MSG_INFO1( spxout << std::endl << "Primal ray (name, id, value):" << std::endl; )
            for( int i = 0; i < SoPlexShell.numColsReal(); ++i )
            {
               if ( isNotZero( ray[i], 0.001 * Real(SoPlexShell.rationalParam(SoPlex2::FEASTOL)) ) )
               {
                  rayobjval += ray[i] * objcoef[i];

                  MSG_INFO1( spxout << colnames[ SoPlexShell.colIdReal(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(17)
                                    << std::setprecision( precision )
                                    << ray[i] << std::endl; )
               }
            }
            MSG_INFO1( spxout << "All other variables have zero value (within " << std::setprecision(1) << 0.001 * SoPlexShell.rationalParam(SoPlex2::FEASTOL) << ")." << std::endl; )
            MSG_INFO1( spxout << "Objective change per unit along primal ray is " << rayobjval << "." << std::endl; )
         }
      }
      break;
   case SPxSolver::INFEASIBLE:
      if ( checkMode )
    MSG_INFO1( spxout << "IEXAMP32 LP is infeasible" << std::endl; )
      else
    MSG_INFO1( spxout << "LP is infeasible" << std::endl; )
      if ( print_solution )
      {
         DVector farkasx(SoPlexShell.numRowsReal());

         if( SoPlexShell.getDualfarkasReal(farkasx) )
         {
            DVector proofvec(SoPlexShell.numColsReal());
            double lhs;
            double rhs;

            lhs = 0.0;
            rhs = 0.0;
            proofvec.clear();
            for( int i = 0; i < SoPlexShell.numRowsReal(); ++i )
            {
               if ( isNotZero( farkasx[i], 0.001 * Real(SoPlexShell.rationalParam(SoPlex2::OPTTOL)) ) )
               {
                  MSG_INFO1( spxout << rownames[ SoPlexShell.rowIdReal(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(16)
                                    << std::setprecision( precision )
                                    << farkasx[i] << "\t"; )
                  LPRow row;
                  SoPlexShell.getRowReal(i, row);
                  if( row.lhs() > -soplex::infinity )
                  {
                     MSG_INFO1( spxout << row.lhs() << " <= "; );
                  }
                  for( int j = 0; j < row.rowVector().size(); ++j )
                  {
                     if( row.rowVector().value(j) > 0 )
                     {
                        MSG_INFO1( spxout << "+"; )
                     }
                     MSG_INFO1( spxout
                        << row.rowVector().value(j) << " "
                        << colnames[ SoPlexShell.colIdReal(row.rowVector().index(j)) ]
                        << " "; );
                  }
                  if( row.rhs() < soplex::infinity )
                  {
                     MSG_INFO1( spxout << "<= " << row.rhs(); );
                  }
                  MSG_INFO1( spxout << std::endl; )
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

            MSG_INFO1( spxout << "All other row multipliers are zero (within " << std::setprecision(1) << 0.001 * SoPlexShell.rationalParam(SoPlex2::OPTTOL) << ")." << std::endl; )
            MSG_INFO1( spxout << "Farkas infeasibility proof: \t"; )
            MSG_INFO1( spxout << lhs << " <= "; )

            bool nonzerofound = false;
            for( int i = 0; i < SoPlexShell.numColsReal(); ++i )
            {
               if ( isNotZero( proofvec[i], 0.001 * Real(SoPlexShell.rationalParam(SoPlex2::OPTTOL)) ) )
               {
                  if( proofvec[i] > 0 )
                  {
                     MSG_INFO1( spxout << "+"; )
                  }
                  MSG_INFO1( spxout << proofvec[i] << " " << colnames[ SoPlexShell.colIdReal(i) ] << " "; )
                  nonzerofound = true;
               }
            }
            if( !nonzerofound )
            {
               MSG_INFO1( spxout << "0 "; );
            }
            MSG_INFO1( spxout << "<= " << rhs << std::endl; );
         }
      }
      if( print_quality )
         displayInfeasibility(SoPlexShell);
      if( write_basis )  // write basis even if we are infeasible
         SoPlexShell.writeBasisFileReal(basisname, &rownames, &colnames);
      break;
   case SPxSolver::ABORT_CYCLING:
      if( checkMode )
    MSG_INFO1( spxout << "EEXAMP40 aborted due to cycling" << std::endl; )
      else
    MSG_INFO1( spxout << "Aborted due to cycling" << std::endl; )
      break;
   case SPxSolver::ABORT_TIME:
      if( checkMode )
    MSG_INFO1( spxout << "IEXAMP33 aborted due to time limit" << std::endl; )
      else
    MSG_INFO1( spxout << "Aborted due to time limit" << std::endl; )
      break;
   case SPxSolver::ABORT_ITER:
      if( checkMode )
    MSG_INFO1( spxout << "IEXAMP34 aborted due to iteration limit" << std::endl; )
      else
    MSG_INFO1( spxout << "Aborted due to iteration limit" << std::endl; )
      break;
   case SPxSolver::ABORT_VALUE:
      if( checkMode )
    MSG_INFO1( spxout << "IEXAMP35 aborted due to objective value limit" << std::endl; )
      else
    MSG_INFO1( spxout << "Aborted due to objective value limit" << std::endl; )
      break;
   case SPxSolver::SINGULAR:
      if( checkMode )
    MSG_INFO1( spxout << "EEXAMP39 basis is singular" << std::endl; )
      else
    MSG_INFO1( spxout << "Basis is singular" << std::endl; )
      break;
   default:
      if( checkMode )
    MSG_INFO1( spxout << "EEXAMP36 An error occurred during " << "the solution process" << std::endl; )
      else
    MSG_INFO1( spxout << "An error occurred during " << "the solution process" << std::endl; )
      break;
   }
   MSG_INFO1( spxout << std::endl; )
}

//------------------------------------------------------------------------
//    main program
//------------------------------------------------------------------------

int main(int argc, char* argv[])
{
   SoPlex2*                  SoPlexShell    = new SoPlex2();
   const char*               filename;
   char*                     basisname      = 0;

   NameSet                   rownames;
   NameSet                   colnames;
   bool                      exact          = false;
   int                       starting       = 0;
   int                       pricing        = 4;
   int                       ratiotest      = 2;
   int                       scaling        = 2;
   int                       simplifying    = 1;
   int                       iterlimit      = -1;
   Real                      timelimit      = -1.0;
   Real                      delta          = DEFAULT_BND_VIOL;
   Real                      feastol        = DEFAULT_BND_VIOL;
   Real                      opttol         = DEFAULT_BND_VIOL;
   Real                      irthreshold    = DEFAULT_BND_VIOL * 1e-6;
   Real                      epsilon        = DEFAULT_EPS_ZERO;
   Real                      epsilon_factor = DEFAULT_EPS_FACTOR;
   Real                      epsilon_update = DEFAULT_EPS_UPDATE;
   int                       verbose        = SPxOut::INFO1;
   bool                      print_solution = false;
   bool                      print_dual     = false;
   bool                      print_quality  = false;
   bool                      read_basis     = false;
   bool                      write_basis    = false;
   bool                      checkMode      = false;
   int                       precision;
   int                       optidx;

   for(optidx = 1; optidx < argc; optidx++)
   {
      if (*argv[optidx] != '-')
         break;

      switch(argv[optidx][1])
      {
         case 'X' :
            exact = true;
            break;
         case 'b' :
            checkParameter(argv[optidx][2], argv); // use -b{r,w}, not -b
            if (argv[optidx][2] == 'r')
               read_basis = true;
            if (argv[optidx][2] == 'w')
               write_basis = true;
            break;
         case 'c' :
            checkParameter(argv[optidx][2], argv); // use -c[0-3], not -c
            setStarter( *SoPlexShell, atoi(&argv[optidx][2]) );
            break;
         case 'd' :
            checkParameter(argv[optidx][2], argv); // use -dx, not -d
            SoPlexShell->setRationalParam( SoPlex2::FEASTOL, atof(&argv[optidx][2]) );
            SoPlexShell->setRationalParam( SoPlex2::OPTTOL, atof(&argv[optidx][2]) );
            break;
         case 'f' :
            checkParameter(argv[optidx][2], argv); // use -fx, not -f
            SoPlexShell->setRationalParam( SoPlex2::FEASTOL, atof(&argv[optidx][2]) );
            break;
         case 'o' :
            checkParameter(argv[optidx][2], argv); // use -ox, not -o
            SoPlexShell->setRationalParam( SoPlex2::OPTTOL, atof(&argv[optidx][2]) );
            break;
         case 'R' :
            checkParameter(argv[optidx][2], argv); // use -Rx, not -R
            SoPlexShell->setRealParam( SoPlex2::IRTHRESHOLD, atof(&argv[optidx][2]) );
            break;
         case 'e':
            SoPlexShell->setIntParam( SoPlex2::ALGORITHM, SoPlex2::ALGORITHM_ENTER );
            break;
         case 'g' :
            checkParameter(argv[optidx][2], argv); // use -g[0-5], not -g
            scaling = atoi(&argv[optidx][2]);
            setScalers( *SoPlexShell, scaling );
            break;
         case 'i' :
            SoPlexShell->setIntParam( SoPlex2::FACTOR_UPDATE_TYPE, SoPlex2::FACTOR_UPDATE_TYPE_ETA );
            break;
         case 'l' :
            checkParameter(argv[optidx][2], argv); // use -lx, not -l
            SoPlexShell->setRealParam( SoPlex2::TIMELIMIT, atoi(&argv[optidx][2]) );
            break;
         case 'L' :
            checkParameter(argv[optidx][2], argv); // use -Lx, not -L
            SoPlexShell->setIntParam( SoPlex2::ITERLIMIT, atoi(&argv[optidx][2]) );
            break;
         case 'p' :
            checkParameter(argv[optidx][2], argv); // use -p[0-5], not -p
            setPricer( *SoPlexShell, atoi(&argv[optidx][2]) );
            break;
         case 'P' :
            SoPlexShell->setBoolParam(SoPlex2::PARTIAL_PRICING, true);
            break;
         case 'q' :
            print_quality = true;
            break;
         case 'r' :
            SoPlexShell->setIntParam( SoPlex2::REPRESENTATION, SoPlex2::REPRESENTATION_ROW );
            break;
         case 's' :
            checkParameter(argv[optidx][2], argv); // use -s[0-1], not -s
            if (argv[optidx][2] >= '0' && argv[optidx][2] <= '9')
               SoPlexShell->setIntParam( SoPlex2::SIMPLIFIER, argv[optidx][2] - '0');
            break;
         case 't' :
            checkParameter(argv[optidx][2], argv); // use -r[0-2], not -r
            setRatiotester( *SoPlexShell, atoi(&argv[optidx][2]) );
            break;
         case 'v' :
            checkParameter(argv[optidx][2], argv); // use -v[0-5], not -v
            if (argv[optidx][2] >= '0' && argv[optidx][2] <= '9')
            {
               verbose = argv[optidx][2] - '0';
               Param::setVerbose( verbose );
            }
            break;
         case 'V' :
            printVersionInfo(checkMode);
            exit(0);
         case 'x' :
            print_solution = true;
            break;
         case 'y' :
            print_dual = true;
            break;
         case 'z' :
            checkParameter(argv[optidx][2], argv); // must not be empty
            checkParameter(argv[optidx][3], argv); // must not be empty
            switch(argv[optidx][2])
            {
            case 'z' :
               SoPlexShell->setRealParam( SoPlex2::EPSILON_ZERO, atof(&argv[optidx][3]) );
               break;
            case 'f' :
               SoPlexShell->setRealParam( SoPlex2::EPSILON_FACTORIZATION, atof(&argv[optidx][3]) );
               break;
            case 'u' :
               SoPlexShell->setRealParam( SoPlex2::EPSILON_UPDATE, atof(&argv[optidx][3]) );
               break;
            default :
               printUsage( argv );
            }
            break;
         case 'C' :
            checkMode = true;
            break;
         case 'h' :
         case '?' :
            printVersionInfo(checkMode);
            //lint -fallthrough
         default :
            printUsage( argv );
      }
   }

   // print version
   printVersionInfo(checkMode);

   // enough arguments?
   if ((argc - optidx) < 1 + (read_basis ? 1 : 0) + (write_basis ? 1 : 0))
      printUsage( argv );
   filename  = argv[optidx];

   ++optidx;

   // switch off simplifier when using a starting basis
   if ( read_basis )
      SoPlexShell->setIntParam( SoPlex2::SIMPLIFIER, SoPlex2::SIMPLIFIER_OFF);

   if ( read_basis || write_basis )
      basisname = strcpy( new char[strlen(argv[optidx]) + 1], argv[optidx] );

   // Set the output precision.
   precision = int(-log10(std::min(feastol, opttol))) + 1;

   std::cout.setf( std::ios::scientific | std::ios::showpoint );
   std::cerr.setf( std::ios::scientific | std::ios::showpoint );

#ifdef SEND_ALL_OUTPUT_TO_FILES
   // Example of redirecting output to different files.
   // Default is cerr for errors and warnings, cout for everything else.
   std::ofstream  myerrstream ( "errwarn.txt" );
   std::ofstream  myinfostream( "infos.txt" );
   redirectOutput(myerrstream, myinfostream);
#endif

   printAlgorithmParameters( *SoPlexShell, checkMode );

   try
   {
      if( exact )
      {
         // read the LP from an input file (.lp or .mps)
         SoPlexShell->readFileRational( filename, &rownames, &colnames );

         // read a basis file if specified
         if (read_basis)
            SoPlexShell->readBasisFileRational( filename, &rownames, &colnames );

         // solve the LP
         SoPlexShell->solveRational();

         // print solution, status, infeasibility system,...
         printSolutionAndStatus( *SoPlexShell, rownames, colnames, precision, print_quality,
                                 print_solution, print_dual, write_basis, basisname, checkMode );
      }
      else
      {
         if ( checkMode )
            MSG_INFO1( spxout << "IEXAMP22 "; )
         MSG_INFO1( spxout << "Loading LP file " << filename << std::endl; )

         // set up timer
         Timer timer;
         timer.start();
         // read the LP from an input file (.lp or .mps)
         if( !SoPlexShell->readFileReal( filename, &rownames, &colnames ) )
         {
            if ( checkMode )
               MSG_INFO1( spxout << "EEXAMP23 "; )
            MSG_INFO1( spxout << "error while reading file \""  << filename << "\"" << std::endl; )
            exit(1);
         }
         timer.stop();

         if ( checkMode )
            MSG_INFO1( spxout << "IEXAMP24 "; )

         MSG_INFO1( spxout << "LP has "
            << SoPlexShell->numRowsReal() << " rows "
            << SoPlexShell->numColsReal() << " columns "
            << SoPlexShell->numNonzerosReal() << " nonzeros"
            << std::endl; )

         if( checkMode )
            MSG_INFO1( spxout << "IEXAMP41 "; )

         MSG_INFO1( std::streamsize prec = spxout.precision();
            spxout << "LP reading time: "
            << std::fixed << std::setprecision(2) << timer.userTime()
            << std::scientific << std::setprecision(int(prec)) << std::endl << std::endl; )

         // read a basis file if specified
         if (read_basis)
         {
            if( !SoPlexShell->readBasisFileReal( filename, &rownames, &colnames ) )
            {
               if ( checkMode )
                  MSG_INFO1( spxout << "EEXAMP25 "; )
               MSG_INFO1( spxout << "error while reading file \""  << filename << "\"" << std::endl; )
               exit(1);
            }
         }

         // solve the LP
         SoPlexShell->solveReal();

         // print solution, status, infeasibility system,...
         printSolutionAndStatus( *SoPlexShell, rownames, colnames, precision, print_quality,
                                 print_solution, print_dual, write_basis, basisname, checkMode );
      }

      // clean up
      delete [] basisname;
      delete SoPlexShell;

      return 0;
   }
   catch(SPxException& x)
   {
      std::cout << "exception caught : " << x.what() << std::endl;
      delete [] basisname;
      delete SoPlexShell;
   }
}
