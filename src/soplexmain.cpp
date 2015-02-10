/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  soplexmain.cpp
 * @brief Command line interface of SoPlex LP solver
 */
#ifndef SOPLEX_LEGACY
#include <assert.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include "soplex.h"
#include "spxgithash.h"
#include "timerfactory.h"

#ifdef SOPLEX_WITH_EGLIB
extern "C" {
#include "EGlib.h"
}
#else
#define EGlpNumStart() {}
#define EGlpNumClear() {}
#endif

using namespace soplex;

// prints usage and command line options
static
void printUsage(const char* const argv[], int idx)
{
   const char* usage =
      "general options:\n"
      "  --readbas=<basfile>    read starting basis from file\n"
      "  --writebas=<basfile>   write terminal basis to file\n"
      "  --writefile=<lpfile>   write LP to file in LP or MPS format depending on extension\n"
      "  --<type>:<name>=<val>  change parameter value using syntax of settings file entries\n"
      "  --loadset=<setfile>    load parameters from settings file (overruled by command line parameters)\n"
      "  --saveset=<setfile>    save parameters to settings file\n"
      "  --diffset=<setfile>    save modified parameters to settings file\n"
      "\n"
      "limits and tolerances:\n"
      "  -t<s>                  set time limit to <s> seconds\n"
      "  -i<n>                  set iteration limit to <n>\n"
      "  -f<eps>                set primal feasibility tolerance to <eps>\n"
      "  -o<eps>                set dual feasibility (optimality) tolerance to <eps>\n"
      "\n"
      "algorithmic settings (* indicates default):\n"
      "  --readmode=<value>     choose reading mode for <lpfile> (0* - floating-point, 1 - rational)\n"
      "  --solvemode=<value>    choose solving mode (0 - floating-point solve, 1* - auto, 2 - force iterative refinement)\n"
      "  -s<value>              choose simplifier/presolver (0 - off, 1* - auto)\n"
      "  -g<value>              choose scaling (0 - off, 1 - uni-equilibrium, 2* - bi-equilibrium, 3 - geometric, 4 - iterated geometric)\n"
      "  -p<value>              choose pricing (0* - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)\n"
      "  -r<value>              choose ratio tester (0 - textbook, 1 - harris, 2 - fast, 3* - boundflipping)\n"
      "\n"
      "display options:\n"
      "  -v<level>              set verbosity to <level> (0 - error, 3 - normal, 5 - high)\n"
      "  -x                     print primal solution\n"
      "  -y                     print dual multipliers\n"
      "  -q                     display detailed statistics\n"
      "  -c                     perform final check of optimal solution in original problem\n"
      "\n";

   if( idx <= 0 )
      std::cerr << "missing input file\n\n";
   else
      std::cerr << "invalid option \"" << argv[idx] << "\"\n\n";

   std::cerr << "usage: " << argv[0] << " " << "[options] <lpfile>\n"
#ifdef SOPLEX_WITH_ZLIB
             << "  <lpfile>               linear program as .mps[.gz] or .lp[.gz] file\n\n"
#else
             << "  <lpfile>               linear program as .mps or .lp file\n\n"
#endif
             << usage;
}

// cleans up C strings
static
void freeStrings(char*& s1, char*& s2, char*& s3, char*& s4, char*& s5)
{
   if( s1 != 0 )
   {
      delete [] s1;
      s1 = 0;
   }
   if( s2 != 0 )
   {
      delete [] s2;
      s2 = 0;
   }
   if( s3 != 0 )
   {
      delete [] s3;
      s3 = 0;
   }
   if( s4 != 0 )
   {
      delete [] s4;
      s4 = 0;
   }
   if( s5 != 0 )
   {
      delete [] s5;
      s5 = 0;
   }
}

/// performs external feasibility check with real type
///@todo implement external check; currently we use the internal methods for convenience
static
void checkSolutionReal(SoPlex& soplex)
{
   if( soplex.hasPrimal() )
   {
      Real boundviol;
      Real rowviol;
      Real sumviol;

      if( soplex.getBoundViolationReal(boundviol, sumviol) && soplex.getRowViolationReal(rowviol, sumviol) )
      {
         MSG_INFO1( soplex.spxout,
            Real maxviol = boundviol > rowviol ? boundviol : rowviol;
            bool feasible = (maxviol <= soplex.realParam(SoPlex::FEASTOL));
            soplex.spxout << "Primal solution " << (feasible ? "feasible" : "infeasible") << " in original problem (max. violation = " << maxviol << ").\n"
            );
      }
      else
      {
         MSG_INFO1( soplex.spxout, soplex.spxout << "Could not check primal solution.\n" );
      }
   }
   else
   {
      MSG_INFO1( soplex.spxout, soplex.spxout << "No primal solution available.\n" );
   }

   if( soplex.hasDual() )
   {
      Real redcostviol;
      Real dualviol;
      Real sumviol;

      if( soplex.getRedCostViolationReal(redcostviol, sumviol) && soplex.getDualViolationReal(dualviol, sumviol) )
      {
         MSG_INFO1( soplex.spxout,
            Real maxviol = redcostviol > dualviol ? redcostviol : dualviol;
            bool feasible = (maxviol <= soplex.realParam(SoPlex::OPTTOL));
            soplex.spxout << "Dual solution " << (feasible ? "feasible" : "infeasible") << " in original problem (max. violation = " << maxviol << ").\n"
            );
      }
      else
      {
         MSG_INFO1( soplex.spxout, soplex.spxout << "Could not check dual solution.\n" );
      }
   }
   else
   {
      MSG_INFO1( soplex.spxout, soplex.spxout << "No dual solution available.\n" );
   }
}

/// performs external feasibility check with rational type
///@todo implement external check; currently we use the internal methods for convenience
static
void checkSolutionRational(SoPlex& soplex)
{
   if( soplex.hasPrimal() )
   {
      Rational boundviol;
      Rational rowviol;
      Rational sumviol;

      if( soplex.getBoundViolationRational(boundviol, sumviol) && soplex.getRowViolationRational(rowviol, sumviol) )
      {
         MSG_INFO1( soplex.spxout,
            Rational maxviol = boundviol > rowviol ? boundviol : rowviol;
            bool feasible = (maxviol <= soplex.realParam(SoPlex::FEASTOL));
            soplex.spxout << "Primal solution " << (feasible ? "feasible" : "infeasible") << " in original problem (max. violation = " << rationalToString(maxviol) << ").\n"
            );
      }
      else
      {
         MSG_INFO1( soplex.spxout, soplex.spxout << "Could not check primal solution.\n" );
      }
   }
   else
   {
      MSG_INFO1( soplex.spxout, soplex.spxout << "No primal solution available.\n" );
   }

   if( soplex.hasDual() )
   {
      Rational redcostviol;
      Rational dualviol;
      Rational sumviol;

      if( soplex.getRedCostViolationRational(redcostviol, sumviol) && soplex.getDualViolationRational(dualviol, sumviol) )
      {
         MSG_INFO1( soplex.spxout,
            Rational maxviol = redcostviol > dualviol ? redcostviol : dualviol;
            bool feasible = (maxviol <= soplex.realParam(SoPlex::OPTTOL));
            soplex.spxout << "Dual solution " << (feasible ? "feasible" : "infeasible") << " in original problem (max. violation = " << rationalToString(maxviol) << ").\n"
            );
      }
      else
      {
         MSG_INFO1( soplex.spxout, soplex.spxout << "Could not check dual solution.\n" );
      }
   }
   else
   {
      MSG_INFO1( soplex.spxout, soplex.spxout << "No dual solution available.\n" );
   }
}

/// performs external feasibility check according to check mode
static
void checkSolution(SoPlex& soplex)
{
   if( soplex.intParam(SoPlex::CHECKMODE) == SoPlex::CHECKMODE_RATIONAL
      || (soplex.intParam(SoPlex::CHECKMODE) == SoPlex::CHECKMODE_AUTO
         && soplex.intParam(SoPlex::READMODE) == SoPlex::READMODE_RATIONAL) )
   {
      checkSolutionRational(soplex);
   }
   else
   {
      checkSolutionReal(soplex);
   }

   MSG_INFO1( soplex.spxout, soplex.spxout << "\n" );
}

/// runs SoPlex command line
int main(int argc, char* argv[])
{
   ///@todo the EGlib version info should be printed after the SoPlex version info
   // initialize EGlib's GMP memory management before any rational numbers are created
   EGlpNumStart();

   SoPlex* soplex;
   NameSet rownames;
   NameSet colnames;
   Timer* readingTime;
   int optidx;

   const char* lpfilename;
   char* readbasname = 0;
   char* writebasname = 0;
   char* writefilename = 0;
   char* loadsetname = 0;
   char* savesetname = 0;
   char* diffsetname = 0;
   bool printPrimal = false;
   bool printDual = false;
   bool displayStatistics = false;
   bool checkSol = false;

   int returnValue = 0;

   // create default timer (CPU time)
   readingTime = TimerFactory::createTimer(Timer::USER_TIME);
   soplex = 0;
   spx_alloc(soplex);
   new (soplex) SoPlex();

   soplex->printVersion();
   MSG_INFO1( soplex->spxout, soplex->spxout << "Copyright (c) 1996-2015 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)\n\n" );

   try
   {
      // no options were given
      if( argc <= 1 )
      {
         printUsage(argv, 0);
         returnValue = 1;
         goto TERMINATE;
      }

      // read arguments from command line
      for( optidx = 1; optidx < argc; optidx++ )
      {
         char* option = argv[optidx];

         // we reached <lpfile>
         if( option[0] != '-' )
            break;

         // option string must start with '-', must contain at least two characters, and exactly two characters if and
         // only if it is -x, -y, -q, or -c
         if( option[0] != '-' || option[1] == '\0'
            || ((option[2] == '\0') != (option[1] == 'x' || option[1] == 'y' || option[1] == 'q' || option[1] == 'c')) )
         {
            printUsage(argv, optidx);
            returnValue = 1;
            goto TERMINATE_FREESTRINGS;
         }

         switch( option[1] )
         {
         case '-' :
            {
               option = &option[2];

               // --readbas=<basfile> : read starting basis from file
               if( strncmp(option, "readbas=", 8) == 0 )
               {
                  if( readbasname == 0 )
                  {
                     char* filename = &option[8];
                     readbasname = strncpy(new char[strlen(filename) + 1], filename, strlen(filename) + 1);
                  }
               }
               // --writebas=<basfile> : write terminal basis to file
               else if( strncmp(option, "writebas=", 9) == 0 )
               {
                  if( writebasname == 0 )
                  {
                     char* filename = &option[9];
                     writebasname = strncpy(new char[strlen(filename) + 1], filename, strlen(filename) + 1);
                  }
               }
               // --writefile=<lpfile> : write LP to file
               else if( strncmp(option, "writefile=", 10) == 0 )
               {
                  if( writefilename == 0 )
                  {
                     char* filename = &option[10];
                     writefilename = strncpy(new char[strlen(filename) + 1], filename, strlen(filename) + 1);
                  }
               }
               // --loadset=<setfile> : load parameters from settings file
               else if( strncmp(option, "loadset=", 8) == 0 )
               {
                  if( loadsetname == 0 )
                  {
                     char* filename = &option[8];
                     loadsetname = strncpy(new char[strlen(filename) + 1], filename, strlen(filename) + 1);
                     if( !soplex->loadSettingsFile(loadsetname) )
                     {
                        printUsage(argv, optidx);
                        returnValue = 1;
                        goto TERMINATE_FREESTRINGS;
                     }
                     else
                     {
                        // we need to start parsing again because some command line parameters might have been overwritten
                        optidx = 0;
                     }
                  }
               }
               // --saveset=<setfile> : save parameters to settings file
               else if( strncmp(option, "saveset=", 8) == 0 )
               {
                  if( savesetname == 0 )
                  {
                     char* filename = &option[8];
                     savesetname = strncpy(new char[strlen(filename) + 1], filename, strlen(filename) + 1);
                  }
               }
               // --diffset=<setfile> : save modified parameters to settings file
               else if( strncmp(option, "diffset=", 8) == 0 )
               {
                  if( diffsetname == 0 )
                  {
                     char* filename = &option[8];
                     diffsetname = strncpy(new char[strlen(filename) + 1], filename, strlen(filename) + 1);
                  }
               }
               // --readmode=<value> : choose reading mode for <lpfile> (0* - floating-point, 1 - rational)
               else if( strncmp(option, "readmode=", 9) == 0 )
               {
                  if( !soplex->setIntParam(SoPlex::READMODE, option[9] - '0') )
                  {
                     printUsage(argv, optidx);
                     returnValue = 1;
                     goto TERMINATE_FREESTRINGS;
                  }
               }
               // --solvemode=<value> : choose solving mode (0* - floating-point solve, 1 - auto, 2 - force iterative refinement)
               else if( strncmp(option, "solvemode=", 10) == 0 )
               {
                  if( !soplex->setIntParam(SoPlex::SOLVEMODE, option[10] - '0') )
                  {
                     printUsage(argv, optidx);
                     returnValue = 1;
                     goto TERMINATE_FREESTRINGS;
                  }
                  // if the LP is parsed rationally and might be solved rationally, we choose automatic syncmode such that
                  // the rational LP is kept after reading
                  else if( soplex->intParam(SoPlex::READMODE) == SoPlex::READMODE_RATIONAL
                     && soplex->intParam(SoPlex::SOLVEMODE) != SoPlex::SOLVEMODE_REAL )
                  {
                     soplex->setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
                  }
               }
               // --<type>:<name>=<val> :  change parameter value using syntax of settings file entries
               else if( !soplex->parseSettingsString(option) )
               {
                  printUsage(argv, optidx);
                  returnValue = 1;
                  goto TERMINATE_FREESTRINGS;
               }
               break;
            }

         case 't' :
            // -t<s> : set time limit to <s> seconds
            if( !soplex->setRealParam(SoPlex::TIMELIMIT, atoi(&option[2])) )
            {
               printUsage(argv, optidx);
               returnValue = 1;
               goto TERMINATE_FREESTRINGS;
            }
            break;

         case 'i' :
            // -i<n> : set iteration limit to <n>
            if( !soplex->setIntParam(SoPlex::ITERLIMIT, atoi(&option[2])) )
            {
               printUsage(argv, optidx);
               returnValue = 1;
               goto TERMINATE_FREESTRINGS;
            }
            break;

         case 'f' :
            // -f<eps> : set primal feasibility tolerance to <eps>
            if( !soplex->setRealParam(SoPlex::FEASTOL, atof(&option[2])) )
            {
               printUsage(argv, optidx);
               returnValue = 1;
               goto TERMINATE_FREESTRINGS;
            }
            break;

         case 'o' :
            // -o<eps> : set dual feasibility (optimality) tolerance to <eps>
            if( !soplex->setRealParam(SoPlex::OPTTOL, atof(&option[2])) )
            {
               printUsage(argv, optidx);
               returnValue = 1;
               goto TERMINATE_FREESTRINGS;
            }
            break;

         case 's' :
            // -s<value> : choose simplifier/presolver (0 - off, 1* - auto)
            if( !soplex->setIntParam(SoPlex::SIMPLIFIER, option[2] - '0') )
            {
               printUsage(argv, optidx);
               returnValue = 1;
               goto TERMINATE_FREESTRINGS;
            }
            break;

         case 'g' :
            // -g<value> : choose scaling (0 - off, 1 - uni-equilibrium, 2* - bi-equilibrium, 3 - geometric, 4 - iterated geometric)
            if( !soplex->setIntParam(SoPlex::SCALER, option[2] - '0') )
            {
               printUsage(argv, optidx);
               returnValue = 1;
               goto TERMINATE_FREESTRINGS;
            }
            break;

         case 'p' :
            // -p<value> : choose pricing (0* - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)
            if( !soplex->setIntParam(SoPlex::PRICER, option[2] - '0') )
            {
               printUsage(argv, optidx);
               returnValue = 1;
               goto TERMINATE_FREESTRINGS;
            }
            break;

         case 'r' :
            // -r<value> : choose ratio tester (0 - textbook, 1 - harris, 2* - fast, 3 - boundflipping)
            if( !soplex->setIntParam(SoPlex::RATIOTESTER, option[2] - '0') )
            {
               printUsage(argv, optidx);
               returnValue = 1;
               goto TERMINATE_FREESTRINGS;
            }
            break;

         case 'v' :
            // -v<level> : set verbosity to <level> (0 - error, 3 - normal, 5 - high)
            if( !soplex->setIntParam(SoPlex::VERBOSITY, option[2] - '0') )
            {
               printUsage(argv, optidx);
               returnValue = 1;
               goto TERMINATE_FREESTRINGS;
            }
            break;

         case 'x' :
            // -x : print primal solution
            printPrimal = true;
            break;

         case 'y' :
            // -y : print dual multipliers
            printDual = true;
            break;

         case 'q' :
            // -q : display detailed statistics
            displayStatistics = true;
            break;

         case 'c' :
            // -c : perform final check of optimal solution in original problem
            checkSol = true;
            break;

            //lint -fallthrough
         default :
            {
               printUsage(argv, optidx);
               returnValue = 1;
               goto TERMINATE_FREESTRINGS;
            }
         }
      }

      MSG_INFO1( soplex->spxout, soplex->printUserSettings(); )

      // no LP file was given and no settings files are written
      if( optidx >= argc && savesetname == 0 && diffsetname == 0 )
      {
         printUsage(argv, 0);
         returnValue = 1;
         goto TERMINATE_FREESTRINGS;
      }

      // ensure that syncmode is not manual
      if( soplex->intParam(SoPlex::SYNCMODE) == SoPlex::SYNCMODE_MANUAL )
      {
         MSG_ERROR( std::cerr << "Error: manual synchronization is invalid on command line.  Change parameter int:syncmode.\n" );
         returnValue = 1;
         goto TERMINATE_FREESTRINGS;
      }

      // save settings files
      if( savesetname != 0 )
      {
         MSG_INFO1( soplex->spxout, soplex->spxout << "Saving parameters to settings file <" << savesetname << "> . . .\n" );
         if( !soplex->saveSettingsFile(savesetname, false) )
         {
            MSG_ERROR( std::cerr << "Error writing parameters to file <" << savesetname << ">\n" );
         }
      }
      if( diffsetname != 0 )
      {
         MSG_INFO1( soplex->spxout, soplex->spxout << "Saving modified parameters to settings file <" << diffsetname << "> . . .\n" );
         if( !soplex->saveSettingsFile(diffsetname, true) )
         {
            MSG_ERROR( std::cerr << "Error writing modified parameters to file <" << diffsetname << ">\n" );
         }
      }

      // no LP file given: exit after saving settings
      if( optidx >= argc )
      {
         if( loadsetname != 0 || savesetname != 0 || diffsetname != 0 )
         {
            MSG_INFO1( soplex->spxout, soplex->spxout << "\n" );
         }
         goto TERMINATE_FREESTRINGS;
      }

      // measure time for reading LP file and basis file
      readingTime->start();

      // if the LP is parsed rationally and might be solved rationally, we choose automatic syncmode such that
      // the rational LP is kept after reading
      if( soplex->intParam(SoPlex::READMODE) == SoPlex::READMODE_RATIONAL
         && soplex->intParam(SoPlex::SOLVEMODE) != SoPlex::SOLVEMODE_REAL )
      {
         soplex->setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
      }

      // read LP from input file
      lpfilename = argv[optidx];
      MSG_INFO1( soplex->spxout, soplex->spxout << "Reading "
         << (soplex->intParam(SoPlex::READMODE) == SoPlex::READMODE_REAL ? "(real)" : "(rational)")
         << " LP file <" << lpfilename << "> . . .\n" );

      if( !soplex->readFile(lpfilename, &rownames, &colnames) )
      {
         MSG_ERROR( std::cerr << "Error while reading file <" << lpfilename << ">.\n" );
         returnValue = 1;
         goto TERMINATE_FREESTRINGS;
      }

      // write LP if specified
      if( writefilename != 0 )
      {
         if( !soplex->writeFileReal(writefilename, &rownames, &colnames) )
         {
            MSG_ERROR( std::cerr << "Error while writing file <" << writefilename << ">.\n\n" );
            returnValue = 1;
            goto TERMINATE_FREESTRINGS;
         }
         else
         {
            MSG_INFO1( soplex->spxout, soplex->spxout << "Written LP to file <" << writefilename << ">.\n\n" );
         }
      }

      // read basis file if specified
      if( readbasname != 0 )
      {
         MSG_INFO1( soplex->spxout, soplex->spxout << "Reading basis file <" << readbasname << "> . . . " );
         if( !soplex->readBasisFile(readbasname, &rownames, &colnames) )
         {
            MSG_ERROR( std::cerr << "Error while reading file <" << readbasname << ">.\n" );
            returnValue = 1;
            goto TERMINATE_FREESTRINGS;
         }
      }

      readingTime->stop();

      MSG_INFO1( soplex->spxout,
         std::streamsize prec = soplex->spxout.precision();
         soplex->spxout << "Reading took "
         << std::fixed << std::setprecision(2) << readingTime->time()
         << std::scientific << std::setprecision(int(prec))
         << " seconds.\n\n" );

      MSG_INFO1( soplex->spxout, soplex->spxout << "LP has " << soplex->numRowsReal() << " rows "
         << soplex->numColsReal() << " columns and " << soplex->numNonzerosReal() << " nonzeros.\n\n" );

      // solve the LP
      soplex->solve();

      // print solution, check solution, and display statistics
      if( printPrimal )
      {
         DVector primal(soplex->numColsReal());
         if( soplex->getPrimalReal(primal) )
         {
            MSG_INFO1( soplex->spxout, soplex->spxout << "\nPrimal solution (name, value):\n"; )
            for( int i = 0; i < soplex->numColsReal(); ++i )
            {
               if ( isNotZero( primal[i] ) )
                  MSG_INFO1( soplex->spxout, soplex->spxout << colnames[i] << "\t"
                                    << std::setw(17)
                                    << std::setprecision(9)
                                    << primal[i] << std::endl; )
            }
            MSG_INFO1( soplex->spxout, soplex->spxout << "All other variables are zero (within "
                              << std::setprecision(1) << std::scientific << Param::epsilon()
                              << std::setprecision(8) << std::fixed << ")." << std::endl; )
         }
         else
            MSG_INFO1( soplex->spxout, soplex->spxout << "No primal solution available.")
      }

      if( printDual )
      {
         DVector dual(soplex->numRowsReal());
         if( soplex->getDualReal(dual) )
         {
            MSG_INFO1( soplex->spxout, soplex->spxout << "\nDual multipliers (name, value):\n"; )
            for( int i = 0; i < soplex->numRowsReal(); ++i )
            {
               if ( isNotZero( dual[i] ) )
                  MSG_INFO1( soplex->spxout, soplex->spxout << rownames[i] << "\t"
                                    << std::setw(17)
                                    << std::setprecision(9)
                                    << dual[i] << std::endl; )
            }
            MSG_INFO1( soplex->spxout, soplex->spxout << "All other dual values are zero (within "
                              << std::setprecision(1) << std::scientific << Param::epsilon()
                              << std::setprecision(8) << std::fixed << ")." << std::endl; )
         }
         else
            MSG_INFO1( soplex->spxout, soplex->spxout << "No dual solution available.")
      }

      if( checkSol )
         checkSolution(*soplex);

      if( displayStatistics )
      {
         MSG_INFO1( soplex->spxout, soplex->spxout << "Statistics\n==========\n\n" );
         soplex->printStatistics(soplex->spxout.getStream(SPxOut::INFO1));
      }

      // write basis file if specified
      if( writebasname != 0 )
      {
         if( !soplex->hasBasis() )
         {
            MSG_WARNING( soplex->spxout, soplex->spxout << "No basis information available.  Could not write file <" << writebasname << ">\n\n" );
         }
         else if( !soplex->writeBasisFile(writebasname, &rownames, &colnames) )
         {
            MSG_ERROR( std::cerr << "Error while writing file <" << writebasname << ">.\n\n" );
            returnValue = 1;
            goto TERMINATE_FREESTRINGS;
         }
         else
         {
            MSG_INFO1( soplex->spxout, soplex->spxout << "Written basis information to file <" << writebasname << ">.\n\n" );
         }
      }
   }
   catch( const SPxException& x )
   {
      MSG_ERROR( std::cerr << "Exception caught: " << x.what() << "\n" );
      returnValue = 1;
      goto TERMINATE_FREESTRINGS;
   }

 TERMINATE_FREESTRINGS:
   freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);

 TERMINATE:
   // because EGlpNumClear() calls mpq_clear() for all mpq_t variables, we need to destroy all objects of class Rational
   // beforehand; hence all Rational objects and all data that uses Rational objects must be allocated dynamically via
   // spx_alloc() and freed here; disabling the list memory is crucial
   soplex->~SoPlex();
   spx_free(soplex);
   Rational::disableListMem();
   EGlpNumClear();
   readingTime->~Timer();
   spx_free(readingTime);

   return returnValue;
}
#else
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
// for simplicity: store whether we are in check mode:
static bool checkMode = false;
//------------------------------------------------------------------------


//------------------------------------------------------------------------
//    class MySoPlex
//------------------------------------------------------------------------

/** LP solver class for the command line. */
class MySoPlex : public SoPlex
{
public:
   /// default constructor
   MySoPlex( SPxOut&                   outstream,
             SPxSolver::Type           p_type = SPxSolver::LEAVE,
             SPxSolver::Representation p_rep  = SPxSolver::COLUMN )
      : SoPlex(outstream, p_type, p_rep)
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

      if ( checkMode )
      {
	 MSG_INFO1( (*spxout), (*spxout) << "IEXAMP05 Violations (max/sum)" << std::endl; )

	 m_solver.qualConstraintViolation(maxviol, sumviol);

	 MSG_INFO1( (*spxout), (*spxout) << "IEXAMP06 Constraints      :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 qualConstraintViolation(maxviol, sumviol);

	 MSG_INFO1( (*spxout), (*spxout) << "IEXAMP07       (unscaled) :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 m_solver.qualBoundViolation(maxviol, sumviol);

	 MSG_INFO1( (*spxout), (*spxout) << "IEXAMP08 Bounds           :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 qualBoundViolation(maxviol, sumviol);

	 MSG_INFO1( (*spxout), (*spxout) << "IEXAMP09       (unscaled) :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 if (!m_vanished)
	 {
	    m_solver.qualSlackViolation(maxviol, sumviol);

	    MSG_INFO1( (*spxout), (*spxout) << "IEXAMP10 Slacks           :"
	       << std::setw(16) << maxviol << "  "
	       << std::setw(16) << sumviol << std::endl; )

	    m_solver.qualRedCostViolation(maxviol, sumviol);

	    MSG_INFO1( (*spxout), (*spxout) << "IEXAMP11 Reduced costs    :"
	       << std::setw(16) << maxviol << "  "
	       << std::setw(16) << sumviol << std::endl; )
#if 0
	    MSG_INFO1( (*spxout), (*spxout) << "IEXAMP12 Proven dual bound:"
	       << std::setw(20)
	       << std::setprecision(20)
	       << m_solver.provedDualbound() << std::endl; )
#endif
	  }
      }
      else
      {
	 MSG_INFO1( (*spxout), (*spxout) << "Violations (max/sum)" << std::endl; )

	 m_solver.qualConstraintViolation(maxviol, sumviol);

	 MSG_INFO1( (*spxout), (*spxout) << "Constraints      :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 qualConstraintViolation(maxviol, sumviol);

	 MSG_INFO1( (*spxout), (*spxout) << "      (unscaled) :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 m_solver.qualBoundViolation(maxviol, sumviol);

	 MSG_INFO1( (*spxout), (*spxout) << "Bounds           :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 qualBoundViolation(maxviol, sumviol);

	 MSG_INFO1( (*spxout), (*spxout) << "      (unscaled) :"
	    << std::setw(16) << maxviol << "  "
	    << std::setw(16) << sumviol << std::endl; )

	 if (!m_vanished)
	 {
	    m_solver.qualSlackViolation(maxviol, sumviol);

	    MSG_INFO1( (*spxout), (*spxout) << "Slacks           :"
	       << std::setw(16) << maxviol << "  "
	       << std::setw(16) << sumviol << std::endl; )

	    m_solver.qualRedCostViolation(maxviol, sumviol);

	    MSG_INFO1( (*spxout), (*spxout) << "Reduced costs    :"
	       << std::setw(16) << maxviol << "  "
	       << std::setw(16) << sumviol << std::endl; )
#if 0
	    MSG_INFO1( (*spxout), (*spxout) << "Proven dual bound:"
	       << std::setw(20)
	       << std::setprecision(20)
	       << m_solver.provedDualbound() << std::endl; )
#endif
	  }
      }
   }
   //------------------------------------------------------------------------
   void displayInfeasibility() const
   {
      assert(m_solver.status() == SPxSolver::INFEASIBLE);

#if 0
      if ( checkMode )
      {
	 if( m_solver.isProvenInfeasible() )
	    MSG_INFO1( (*spxout), (*spxout) << "IEXAMP13 Infeasibility is proven." << std::endl; )
         else
	    MSG_INFO1( (*spxout), (*spxout) << "IEXAMP13 Infeasibility could not be proven!" << std::endl; )
      }
      else
      {
         if ( m_solver.isProvenInfeasible() )
	 {
	    MSG_INFO1( (*spxout), (*spxout) << "Infeasibility is proven." << std::endl; )
	 }
	 else
	 {
	    MSG_INFO1( (*spxout), (*spxout) << "Infeasibility could not be proven!" << std::endl; )
	 }
      }
#endif
   }
};


//------------------------------------------------------------------------
//    Helpers
//------------------------------------------------------------------------

static
void print_version_info()
{
   const char* banner1 =
   "************************************************************************\n"
   "*                                                                      *\n"
   "*       SoPlex --- the Sequential object-oriented simPlex.             *\n"
   ;

   const char* banner2 =
   "*                                                                      *\n"
   "*    Copyright (C) 1996-2015 Konrad-Zuse-Zentrum                       *\n"
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
   "* Copyright (C) 1996-2015 Konrad-Zuse-Zentrum                          *\n"
   "*                         fuer Informationstechnik Berlin              *\n"
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
void print_usage_and_exit( const char* const argv[] )
{
   const char* usage =
      "[options] LPfile [Basfile]\n\n"
      "          LPfile can be either in MPS or LPF format\n\n"
      "options:  (*) indicates default\n"
      "          (!) indicates experimental features which may give wrong results\n"
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
      " -zz       set general zero tolerance\n"
      " -zf       set factorization zero tolerance\n"
      " -zu       set update zero tolerance\n"
      " -v        set verbosity Level: from 0 (ERROR) to 5 (INFO3), default 3 (INFO1)\n"
      " -V        show program version\n"
      " -C        check mode (for check scripts)\n"
      " -h        show this help\n\n"
      "Simplifier:  Scaler:           Starter:    Pricer:        Ratiotester:\n"
      " -s0 none     -g0 none          -c0 none*   -p0 Textbook   -t0 Textbook\n"
      " -s1 Main*    -g1 uni-Equi      -c1 Weight  -p1 ParMult    -t1 Harris\n"
      "              -g2 bi-Equi*      -c2 Sum     -p2 Devex      -t2 Fast\n"
      "              -g3 bi-Equi+Geo1  -c3 Vector  -p3 Hybrid!    -t3 Bound Flipping*\n"
      "              -g4 bi-Equi+Geo8              -p4 Steep*\n"
      "                                            -p5 Weight\n"
      "                                            -p6 SteepExactSetup\n"
      ;

   std::cerr << "usage: " << argv[0] << " " << usage << std::endl;
   exit(0);
}

//------------------------------------------------------------------------
static
void check_parameter(const char param, const char* const argv[])
{
   if (param == '\0')
      print_usage_and_exit( argv );
}

//------------------------------------------------------------------------
static
void print_algorithm_parameters(
   MySoPlex&                       work,
   const SPxSolver::Representation representation,
   const SLUFactor::UpdateType     update
   )
{
   if ( checkMode )
   {
      MSG_INFO1( (*work.spxout), (*work.spxout)
	 << "IEXAMP12 Feastol        = "
	 << std::setw(16) << work.feastol() << std::endl
	 << "IEXAMP52 Opttol         = "
	 << std::setw(16) << work.opttol() << std::endl
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
   else
   {
      MSG_INFO1( (*work.spxout), (*work.spxout)
	 << "SoPlex parameters: " << std::endl
	 << "Feastol        = "
	 << std::setw(16) << work.feastol() << std::endl
	 << "Opttol         = "
	 << std::setw(16) << work.opttol() << std::endl
	 << "Epsilon Zero   = "
	 << std::setw(16) << Param::epsilon() << std::endl
	 << "Epsilon Factor = "
	 << std::setw(16) << Param::epsilonFactorization() << std::endl
	 << "Epsilon Update = "
	 << std::setw(16) << Param::epsilonUpdate() << std::endl
	 << std::endl
	 << "algorithm      = " << (work.type() == SPxSolver::ENTER ? "Entering" : "Leaving")
	 << std::endl
	 << "representation = " << (representation == SPxSolver::ROW ? "Row" : "Column")
	 << std::endl
	 << "update         = " << (update == SLUFactor::ETA ? "Eta" : "Forest-Tomlin")
	 << std::endl; )
   }
}

//------------------------------------------------------------------------
static
SPxPricer* get_pricer(const int pricing, SPxOut* spxout)
{
   SPxPricer* pricer = 0;
   switch(pricing)
   {
   case 6 :
      pricer = new SPxSteepExPR;
      break;
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

   assert(pricer != 0);
   if ( checkMode )
#ifdef PARTIAL_PRICING
      MSG_INFO1( (*spxout), (*spxout) << "IEXAMP17 " << pricer->getName() << " pricing"
                        << " (partial, size = " << MAX_PRICING_CANDIDATES << ")"
                        << std::endl; )
#else
      MSG_INFO1( (*spxout), (*spxout) << "IEXAMP17 " << pricer->getName() << " pricing"
                        << std::endl; )
#endif
   else
#ifdef PARTIAL_PRICING
      MSG_INFO1( (*spxout), (*spxout) << "pricing        = " << pricer->getName()
                        << " (partial, size = " << MAX_PRICING_CANDIDATES << ")"
                        << std::endl; )
#else
      MSG_INFO1( (*spxout), (*spxout) << "pricing        = " << pricer->getName()
                        << std::endl; )
#endif
   return pricer;
}

//------------------------------------------------------------------------
static
SPxRatioTester* get_ratio_tester(const int ratiotest, SPxOut* spxout)
{
   SPxRatioTester* ratiotester = 0;
   switch(ratiotest)
   {
   case 3 :
      ratiotester = new SPxBoundFlippingRT;
      break;
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

   assert(ratiotester != 0);
   if ( checkMode )
      MSG_INFO1( (*spxout), (*spxout) << "IEXAMP18 " << ratiotester->getName() << " ratiotest" << std::endl; )
   else
      MSG_INFO1( (*spxout), (*spxout) << "ratiotest      = " << ratiotester->getName() << std::endl; )
   return ratiotester;
}

//------------------------------------------------------------------------
static
void get_scalers(
   SPxScaler*& prescaler,
   SPxScaler*& postscaler,
   const int   scaling,
   SPxOut*     spxout
   )
{
   switch(scaling)
   {
   case 4:
      prescaler  = new SPxEquiliSC(true);
      postscaler = new SPxGeometSC(8);
      break;
   case 3:
      prescaler  = new SPxEquiliSC(true);
      postscaler = new SPxGeometSC(1);
      break;
   case 2 :
      prescaler  = new SPxEquiliSC(true);
      postscaler = 0;
      break;
   case 1 :
      prescaler  = new SPxEquiliSC(false);
      postscaler = 0;
      break;
   case 0 :
      /*FALLTHROUGH*/
   default :
      prescaler  = 0;
      postscaler = 0;
      break;
   }

   if ( checkMode )
   {
      MSG_INFO1( (*spxout), (*spxout) << "IEXAMP19 "
	 << ((prescaler != 0) ? prescaler->getName() : "no")
	 << " / "
	 << ((postscaler != 0) ? postscaler->getName() : "no")
	 << " scaling" << std::endl; )
   }
   else
   {
      MSG_INFO1( (*spxout), (*spxout) << "scaling        = "
	 << ((prescaler != 0) ? prescaler->getName() : "no")
	 << " / "
	 << ((postscaler != 0) ? postscaler->getName() : "no")
	 << std::endl; )
   }
}

//------------------------------------------------------------------------
static
SPxSimplifier* get_simplifier(const int simplifying, SPxOut* spxout)
{
   SPxSimplifier* simplifier = 0;
   switch(simplifying)
   {
   case 1 :
      simplifier = new SPxMainSM;
      break;
   case 0  :
      /*FALLTHROUGH*/
   default :
      assert(simplifier == 0);
      break;
   }

   if ( checkMode )
      MSG_INFO1( (*spxout), (*spxout) << "IEXAMP20 " << ((simplifier == 0) ? "no" : simplifier->getName()) << " simplifier" << std::endl; )
   else
      MSG_INFO1( (*spxout), (*spxout) << "simplifier     = " << ((simplifier == 0) ? "no" : simplifier->getName()) << std::endl; )
   return simplifier;
}

//------------------------------------------------------------------------
static
SPxStarter* get_starter(const int starting, SPxOut* spxout)
{
   SPxStarter* starter = 0;
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

   if ( checkMode )
      MSG_INFO1( (*spxout), (*spxout) << "IEXAMP21 " << ((starter == 0) ? "no" : starter->getName()) << " starter" << std::endl; )
   else
      MSG_INFO1( (*spxout), (*spxout) << "starter        = " << ((starter == 0) ? "no" : starter->getName()) << std::endl; )

   return starter;
}

//------------------------------------------------------------------------
#ifdef SEND_ALL_OUTPUT_TO_FILES
static
void redirect_output(
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
void read_input_file(
   MySoPlex&      work,
   const char*    filename,
   NameSet&       rownames,
   NameSet&       colnames)
{
   if ( checkMode )
      MSG_INFO1( (*work.spxout), (*work.spxout) << "IEXAMP22 loading LP file " << filename << std::endl; )
   else
      MSG_INFO1( (*work.spxout), (*work.spxout) << "\nLoading LP file " << filename << std::endl; )

   UserTimer timer;
   timer.start();

   if ( ! work.readFile(filename, &rownames, &colnames, 0) )
   {
      if ( checkMode )
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "EEXAMP23 error while reading file \"" << filename << "\"" << std::endl; )
      else
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "error while reading file \""  << filename << "\"" << std::endl; )
      exit(1);
   }
   assert(work.isConsistent());

   timer.stop();

   if ( checkMode )
   {
      MSG_INFO1( (*work.spxout), (*work.spxout) << "IEXAMP24 LP has "
	 << work.nRows() << " rows "
	 << work.nCols() << " columns "
	 << work.nNzos() << " nonzeros"
	 << std::endl; )

      MSG_INFO1( (*work.spxout), (*work.spxout) << "IEXAMP41 LP reading time: " << timer.time() << std::endl; )
   }
   else
   {
      MSG_INFO1( (*work.spxout), (*work.spxout) << "LP has "
	 << work.nRows() << " rows "
	 << work.nCols() << " columns "
	 << work.nNzos() << " nonzeros"
	 << std::endl; )

      MSG_INFO1( (*work.spxout),
	 std::streamsize prec = (*work.spxout).precision();
	 (*work.spxout) << "LP reading time: " << std::fixed << std::setprecision(2) << timer.time();
	 (*work.spxout) << std::scientific << std::setprecision(int(prec)) << std::endl; )
   }
}

//------------------------------------------------------------------------
static
void read_basis_file(
   MySoPlex&      work    ,
   const char*    filename,
   const NameSet* rownames,
   const NameSet* colnames)
{
   MSG_INFO1( (*work.spxout), (*work.spxout) << "Reading basis from file (disables simplifier)" << std::endl; )
   if (!work.readBasisFile(filename, rownames, colnames))
   {
      if ( checkMode )
         MSG_INFO1( (*work.spxout), (*work.spxout) << "EEXAMP25 error while reading file \"" << filename << "\"" << std::endl; )
      else
         MSG_INFO1( (*work.spxout), (*work.spxout) << "Error while reading file \"" << filename << "\"" << std::endl; )
      exit(1);
   }
}

//------------------------------------------------------------------------
static
void solve_LP(MySoPlex& work)
{
   UserTimer timer;
   timer.start();

   if ( checkMode )
      MSG_INFO1( (*work.spxout), (*work.spxout) << "IEXAMP26 solving LP" << std::endl; )
   else
      MSG_INFO1( (*work.spxout), (*work.spxout) << "\nSolving LP ..." << std::endl; )

   work.solve();
   timer.stop();

   MSG_INFO1( (*work.spxout), (*work.spxout) << "\nSoPlex statistics:\n" << work.statistics(); )
}

//------------------------------------------------------------------------
static
void print_solution_and_status(
   MySoPlex&            work,
   const NameSet&       rownames,
   const NameSet&       colnames,
   const int            precision,
   const bool           print_quality,
   const bool           print_solution,
   const bool           print_dual,
   const bool           write_basis,
   const char*          basisname
   )
{
   // get the solution status
   SPxSolver::Status stat = work.status();

   if ( ! checkMode )
      MSG_INFO1( (*work.spxout), (*work.spxout) << std::endl; )
   switch (stat)
   {
   case SPxSolver::OPTIMAL:
      if ( checkMode )
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "IEXAMP29 solution value is: " << std::setprecision( precision ) << work.objValue() << std::endl; )
      else
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "Solution value is: " << std::setprecision( precision ) << work.objValue() << std::endl; )

      if ( print_quality )
         work.displayQuality();

      if ( print_solution )
      {
         DVector objx(work.nCols());

         if( work.getPrimal(objx) != SPxSolver::ERROR )
         {
            MSG_INFO1( (*work.spxout), (*work.spxout) << std::endl << "Primal solution (name, id, value):" << std::endl; )
            for( int i = 0; i < work.nCols(); ++i )
            {
               if ( isNotZero( objx[i], 0.001 * work.feastol() ) )
                  MSG_INFO1( (*work.spxout), (*work.spxout) << colnames[ work.cId(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(17)
                                    << std::setprecision( precision )
                                    << objx[i] << std::endl; )
            }
            MSG_INFO1( (*work.spxout), (*work.spxout) << "All other variables are zero (within " << std::setprecision(1) << 0.001*work.feastol() << ")." << std::endl; )
         }
      }
      if ( print_dual )
      {
         DVector objy(work.nRows());
         bool allzero = true;

         if( work.getDual(objy) != SPxSolver::ERROR )
         {
            MSG_INFO1( (*work.spxout), (*work.spxout) << std::endl << "Dual multipliers (name, id, value):" << std::endl; )
            for( int i = 0; i < work.nRows(); ++i )
            {
               if ( isNotZero( objy[i] , 0.001 * work.opttol() ) )
               {
                  MSG_INFO1( (*work.spxout), (*work.spxout) << rownames[ work.rId(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(17)
                                    << std::setprecision( precision )
                                    << objy[i] << std::endl; )
                  allzero = false;
               }
            }

            MSG_INFO1( (*work.spxout), (*work.spxout) << "All " << (allzero ? "" : "other ") << "dual values are zero (within "
                              << std::setprecision(1) << 0.001*work.opttol() << ")." << std::endl; )

            if( !allzero )
            {
               if( work.spxSense() == SPxLP::MINIMIZE )
               {
                  MSG_INFO1( (*work.spxout), (*work.spxout) << "Minimizing: a positive/negative value corresponds to left-hand (>=) resp. right-hand (<=) side."
                                    << std::endl; )
               }
               else
               {
                  MSG_INFO1( (*work.spxout), (*work.spxout) << "Maximizing: a positive/negative value corresponds to right-hand (<=) resp. left-hand (>=) side."
                                    << std::endl; )
               }
            }
         }
      }
      if ( write_basis )
      {
         MSG_INFO1( (*work.spxout), (*work.spxout) << "Writing basis of original problem to file " << basisname << std::endl; )
         if ( ! work.writeBasisFile( basisname, &rownames, &colnames ) )
         {
            if ( checkMode )
               MSG_INFO1( (*work.spxout), (*work.spxout) << "EEXAMP30 error while writing file \"" << basisname << "\"" << std::endl; )
            else
               MSG_INFO1( (*work.spxout), (*work.spxout) << "Error while writing file \"" << basisname << "\"" << std::endl; )
         }
      }
      break;
   case SPxSolver::UNBOUNDED:
      if ( checkMode )
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "IEXAMP31 LP is unbounded" << std::endl; )
      else
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "LP is unbounded" << std::endl; )

      if ( print_solution )
      {
         DVector objx(work.nCols());
         if( work.getPrimal(objx) != SPxSolver::ERROR )
         {
            MSG_INFO1( (*work.spxout), (*work.spxout) << std::endl << "Primal solution (name, id, value):" << std::endl; )
            for( int i = 0; i < work.nCols(); ++i )
            {
               if ( isNotZero( objx[i], 0.001 * work.feastol() ) )
                  MSG_INFO1( (*work.spxout), (*work.spxout) << colnames[ work.cId(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(17)
                                    << std::setprecision( precision )
                                    << objx[i] << std::endl; )
            }
            MSG_INFO1( (*work.spxout), (*work.spxout) << "All other variables are zero (within " << std::setprecision(1) << 0.001*work.feastol() << ")." << std::endl; )
         }

         DVector objcoef(work.nCols());
         DVector ray(work.nCols());
         if( work.getPrimalray(ray) != SPxSolver::ERROR )
         {
            Real rayobjval = 0.0;

            work.getObj(objcoef);

            MSG_INFO1( (*work.spxout), (*work.spxout) << std::endl << "Primal ray (name, id, value):" << std::endl; )
            for( int i = 0; i < work.nCols(); ++i )
            {
               if ( isNotZero( ray[i], 0.001 * work.feastol() ) )
               {
                  rayobjval += ray[i] * objcoef[i];

                  MSG_INFO1( (*work.spxout), (*work.spxout) << colnames[ work.cId(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(17)
                                    << std::setprecision( precision )
                                    << ray[i] << std::endl; )
               }
            }
            MSG_INFO1( (*work.spxout), (*work.spxout) << "All other variables have zero value (within " << std::setprecision(1) << 0.001*work.feastol() << ")." << std::endl; )
            MSG_INFO1( (*work.spxout), (*work.spxout) << "Objective change per unit along primal ray is " << rayobjval << "." << std::endl; )
         }
      }
      break;
   case SPxSolver::INFEASIBLE:
      if ( checkMode )
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "IEXAMP32 LP is infeasible" << std::endl; )
      else
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "LP is infeasible" << std::endl; )
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
               if ( isNotZero( farkasx[i], 0.001 * work.opttol() ) )
               {
                  MSG_INFO1( (*work.spxout), (*work.spxout) << rownames[ work.rId(i) ] << "\t"
                                    << i << "\t"
                                    << std::setw(16)
                                    << std::setprecision( precision )
                                    << farkasx[i] << "\t"; )
                  LPRow row;
                  work.getRow(i, row);
                  if( row.lhs() > -soplex::infinity )
                  {
                     MSG_INFO1( (*work.spxout), (*work.spxout) << row.lhs() << " <= "; );
                  }
                  for( int j = 0; j < row.rowVector().size(); ++j )
                  {
                     if( row.rowVector().value(j) > 0 )
                     {
                        MSG_INFO1( (*work.spxout), (*work.spxout) << "+"; )
                     }
                     MSG_INFO1( (*work.spxout), (*work.spxout)
                        << row.rowVector().value(j) << " "
                        << colnames[ work.cId(row.rowVector().index(j)) ]
                        << " "; );
                  }
                  if( row.rhs() < soplex::infinity )
                  {
                     MSG_INFO1( (*work.spxout), (*work.spxout) << "<= " << row.rhs(); );
                  }
                  MSG_INFO1( (*work.spxout), (*work.spxout) << std::endl; )
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

            MSG_INFO1( (*work.spxout), (*work.spxout) << "All other row multipliers are zero (within " << std::setprecision(1) << 0.001*work.opttol() << ")." << std::endl; )
            MSG_INFO1( (*work.spxout), (*work.spxout) << "Farkas infeasibility proof: \t"; )
            MSG_INFO1( (*work.spxout), (*work.spxout) << lhs << " <= "; )

            bool nonzerofound = false;
            for( int i = 0; i < work.nCols(); ++i )
            {
               if ( isNotZero( proofvec[i], 0.001 * work.opttol() ) )
               {
                  if( proofvec[i] > 0 )
                  {
                     MSG_INFO1( (*work.spxout), (*work.spxout) << "+"; )
                  }
                  MSG_INFO1( (*work.spxout), (*work.spxout) << proofvec[i] << " " << colnames[ work.cId(i) ] << " "; )
                  nonzerofound = true;
               }
            }
            if( !nonzerofound )
            {
               MSG_INFO1( (*work.spxout), (*work.spxout) << "0 "; );
            }
            MSG_INFO1( (*work.spxout), (*work.spxout) << "<= " << rhs << std::endl; );
         }
      }
      if ( print_quality )
         work.displayInfeasibility();
      if ( write_basis )  // write basis even if we are infeasible
         if ( ! work.writeBasisFile( basisname, &rownames, &colnames ) )
         {
	    if ( checkMode )
	       MSG_INFO1( (*work.spxout), (*work.spxout) << "EEXAMP30 error while writing file \"" << basisname << "\"" << std::endl; )
	    else
	       MSG_INFO1( (*work.spxout), (*work.spxout) << "Error while writing file \"" << basisname << "\"" << std::endl; )
         }
      break;
   case SPxSolver::ABORT_CYCLING:
      if ( checkMode )
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "EEXAMP40 aborted due to cycling" << std::endl; )
      else
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "Aborted due to cycling" << std::endl; )
      break;
   case SPxSolver::ABORT_TIME:
      if ( checkMode )
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "IEXAMP33 aborted due to time limit" << std::endl; )
      else
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "Aborted due to time limit" << std::endl; )
      break;
   case SPxSolver::ABORT_ITER:
      if ( checkMode )
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "IEXAMP34 aborted due to iteration limit" << std::endl; )
      else
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "Aborted due to iteration limit" << std::endl; )
      break;
   case SPxSolver::ABORT_VALUE:
      if ( checkMode )
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "IEXAMP35 aborted due to objective value limit" << std::endl; )
      else
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "Aborted due to objective value limit" << std::endl; )
      break;
   case SPxSolver::SINGULAR:
      if ( checkMode )
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "EEXAMP39 basis is singular" << std::endl; )
      else
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "Basis is singular" << std::endl; )
      break;
   default:
      if ( checkMode )
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "EEXAMP36 An error occurred during " << "the solution process" << std::endl; )
      else
	 MSG_INFO1( (*work.spxout), (*work.spxout) << "An error occurred during " << "the solution process" << std::endl; )
      break;
   }
   MSG_INFO1( (*work.spxout), (*work.spxout) << std::endl; )
}

//------------------------------------------------------------------------
static
void clean_up(
   SPxScaler*&       prescaler,
   SPxScaler*&       postscaler,
   SPxSimplifier*&   simplifier,
   SPxStarter*&      starter,
   SPxPricer*&       pricer,
   SPxRatioTester*&  ratiotester,
   char*&            basisname
   )
{
   if ( prescaler != 0 )
   {
      delete prescaler;
      prescaler = 0;
   }
   if ( postscaler != 0 )
   {
      delete postscaler;
      postscaler = 0;
   }
   if ( simplifier != 0 )
   {
      delete simplifier;
      simplifier = 0;
   }
   if ( starter != 0 )
   {
      delete starter;
      starter = 0;
   }

   assert( pricer != 0 );
   delete pricer;
   pricer = 0;

   assert( ratiotester != 0 );
   delete ratiotester;
   ratiotester = 0;

   if ( basisname != 0 )
      delete [] basisname;
   basisname = 0;
}

//------------------------------------------------------------------------
//    main program
//------------------------------------------------------------------------

int main(int argc, char* argv[])
{
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

   SPxOut                    spxout;

   try {
      NameSet                   rownames;
      NameSet                   colnames;
      int                       starting       = 0;
      int                       pricing        = 4;
      int                       ratiotest      = 3;
      int                       scaling        = 2;
      int                       simplifying    = 1;
      int                       iterlimit      = -1;
      Real                      timelimit      = -1.0;
      Real                      delta          = DEFAULT_BND_VIOL;
      Real                      feastol        = DEFAULT_BND_VIOL;
      Real                      opttol         = DEFAULT_BND_VIOL;
      Real                      epsilon        = DEFAULT_EPS_ZERO;
      Real                      epsilon_factor = DEFAULT_EPS_FACTOR;
      Real                      epsilon_update = DEFAULT_EPS_UPDATE;
      int                       verbose        = SPxOut::INFO1;
      bool                      print_solution = false;
      bool                      print_dual     = false;
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
            check_parameter(argv[optidx][2], argv); // use -b{r,w}, not -b
            if (argv[optidx][2] == 'r')
               read_basis = true;
            if (argv[optidx][2] == 'w')
               write_basis = true;
            break;
         case 'c' :
            check_parameter(argv[optidx][2], argv); // use -c[0-3], not -c
            starting = atoi(&argv[optidx][2]);
            break;
         case 'd' :
            check_parameter(argv[optidx][2], argv); // use -dx, not -d
            delta = atof(&argv[optidx][2]);
            break;
         case 'f' :
            check_parameter(argv[optidx][2], argv); // use -fx, not -f
            feastol = atof(&argv[optidx][2]);
            break;
         case 'o' :
            check_parameter(argv[optidx][2], argv); // use -ox, not -o
            opttol = atof(&argv[optidx][2]);
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
         case 'L' :
            if (argv[optidx][2] == '\0' )  // use -Lx, not -L
               print_usage_and_exit( argv );
            iterlimit = atoi(&argv[optidx][2]);
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
            print_version_info();
            exit(0);
         case 'x' :
            print_solution = true;
            break;
         case 'y' :
            print_dual = true;
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
         case 'C' :
            checkMode = true;
            break;
         case 'h' :
         case '?' :
            print_version_info();
            //lint -fallthrough
         default :
            print_usage_and_exit( argv );
         }
      }

      // print version
      print_version_info();

      // enough arguments?
      if ((argc - optidx) < 1 + (read_basis ? 1 : 0) + (write_basis ? 1 : 0))
         print_usage_and_exit( argv );
      filename  = argv[optidx];

      ++optidx;

      // switch off simplifier when using a starting basis
      if ( read_basis )
         simplifying = 0;

      if ( read_basis || write_basis )
         basisname = strcpy( new char[strlen(argv[optidx]) + 1], argv[optidx] );

      // set some algorithm parameters
      Param::setEpsilon             ( epsilon );
      Param::setEpsilonFactorization( epsilon_factor );
      Param::setEpsilonUpdate       ( epsilon_update );
      spxout.setVerbosity           ( (SPxOut::Verbosity) verbose );

      // Set the output precision.
      precision = int(-log10(std::min(feastol, opttol))) + 1;

      std::cout.setf( std::ios::scientific | std::ios::showpoint );
      std::cerr.setf( std::ios::scientific | std::ios::showpoint );

#ifdef SEND_ALL_OUTPUT_TO_FILES
      // Example of redirecting output to different files.
      // Default is cerr for errors and warnings, cout for everything else.
      std::ofstream  myerrstream ( "errwarn.txt" );
      std::ofstream  myinfostream( "infos.txt" );
      redirect_output(myerrstream, myinfostream);
#endif

      // create an instance of MySoPlex
      MySoPlex work( spxout, type, representation );
      work.setOutstream         ( spxout );
      work.setUtype             ( update );
      work.setFeastol           ( std::min(feastol, delta) );
      work.setOpttol            ( std::min(opttol, delta) );
      work.setTerminationTime   ( timelimit );
      work.setTerminationIter   ( iterlimit );
      print_algorithm_parameters( work, representation, update );
      assert( work.isConsistent() );

      // set pricer, starter, simplifier, and ratio tester
      work.setPricer    ( pricer      = get_pricer      (pricing, work.spxout) );
      work.setStarter   ( starter     = get_starter     (starting, work.spxout) );
      work.setSimplifier( simplifier  = get_simplifier  (simplifying, work.spxout) );
      work.setTester    ( ratiotester = get_ratio_tester(ratiotest, work.spxout) );
      assert(work.isConsistent());

      // set pre- and postscaler
      get_scalers(prescaler, postscaler, scaling, work.spxout);
      work.setPreScaler (prescaler);
      work.setPostScaler(postscaler);
      assert(work.isConsistent());

      // read the LP from an input file (.lp or .mps)
      read_input_file(work, filename, rownames, colnames);

      // read a basis file if specified
      if (read_basis)
         read_basis_file(work, basisname, &rownames, &colnames);

      // solve the LP
      solve_LP(work);

      // print solution, status, infeasibility system,...
      print_solution_and_status(work, rownames, colnames, precision, print_quality,
                                print_solution, print_dual, write_basis, basisname);

      // clean up
      clean_up(prescaler, postscaler, simplifier, starter, pricer, ratiotester, basisname);

      return 0;
   }
   catch( const SPxException& x )
   {
      std::cout << "exception caught : " << x.what() << std::endl;
      delete [] basisname;
      if (simplifier)
         delete simplifier;
      delete starter;
      delete pricer;
      delete ratiotester;
      delete prescaler;
      delete postscaler;
   }
}
#endif
