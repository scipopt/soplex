/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2014 Konrad-Zuse-Zentrum                            */
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

#include <assert.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include "soplex.h"
#include "spxgithash.h"
#include "timer.h"

using namespace soplex;


/// prints version and compilation options
static
void printVersionInfo()
{
   MSG_INFO1( spxout << "SoPlex version " << SOPLEX_VERSION/100
      << "." << (SOPLEX_VERSION % 100)/10
      << "." << SOPLEX_VERSION % 10
#if SOPLEX_SUBVERSION > 0
      << "." << SOPLEX_SUBVERSION
#endif
#ifndef NDEBUG
      << " [mode: debug]"
#else
      << " [mode: optimized]"
#endif
      << " [precision: " << (int)sizeof(Real) << " byte]"
#ifdef SOPLEX_WITH_GMP
      << " [rational: gmp]"
#elif SOPLEX_WITH_GMPXX
      << " [rational: gmpxx]"
#else
      << " [rational: long double]"
#endif
      << " [githash: " << std::setw(13) << std::setiosflags(std::ios::left) << getGitHash() << "]\n" );

   MSG_INFO1( spxout << "Copyright (c) 1996-2014 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)\n\n" );
}

// prints usage and command line options
static
void printUsage(const char* const argv[], int idx)
{
   const char* usage =
      "general options:\n"
      "  --readbas=<basfile>    read starting basis from file\n"
      "  --writebas=<basfile>   write terminal basis to file\n"
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
      "  --solvemode=<value>    choose solving mode (0* - floating-point solve, 1 - auto, 2 - force iterative refinement)\n"
      "  -s<value>              choose simplifier/presolver (0 - off, 1* - auto)\n"
      "  -g<value>              choose scaling (0 - off, 1 - uni-equilibrium, 2* - bi-equilibrium, 3 - geometric, 4 - iterated geometric)\n"
      "  -p<value>              choose pricing (0* - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)\n"
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
#ifdef WITH_ZLIB
             << "  <lpfile>               linear program as .mps[.gz] or .lp[.gz] file\n\n"
#else
             << "  <lpfile>               linear program as .mps or .lp file\n\n"
#endif
             << usage;

   exit(1);
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

/// runs SoPlex command line
int main(int argc, char* argv[])
{
   SoPlex soplex;
   NameSet rownames;
   NameSet colnames;
   Timer readingTimer;
   int optidx;

   const char* lpfilename;
   char* readbasname = 0;
   char* writebasname = 0;
   char* loadsetname = 0;
   char* savesetname = 0;
   char* diffsetname = 0;
   bool printPrimal = false;
   bool printDual = false;
   bool displayStatistics = false;
   bool checkSol = false;

   printVersionInfo();

   try
   {
      // no options were given
      if( argc <= 1 )
         printUsage(argv, 0);

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
            freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
            printUsage(argv, optidx);
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
                     readbasname = strcpy(new char[strlen(filename) + 1], filename);
                  }
               }
               // --writebas=<basfile> : write terminal basis to file
               else if( strncmp(option, "writebas=", 9) == 0 )
               {
                  if( writebasname == 0 )
                  {
                     char* filename = &option[9];
                     writebasname = strcpy(new char[strlen(filename) + 1], filename);
                  }
               }
               // --loadset=<setfile> : load parameters from settings file
               else if( strncmp(option, "loadset=", 8) == 0 )
               {
                  if( loadsetname == 0 )
                  {
                     char* filename = &option[8];
                     loadsetname = strcpy(new char[strlen(filename) + 1], filename);
                     if( !soplex.loadSettingsFile(loadsetname) )
                     {
                        freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
                        printUsage(argv, optidx);
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
                     savesetname = strcpy(new char[strlen(filename) + 1], filename);
                  }
               }
               // --diffset=<setfile> : save modified parameters to settings file
               else if( strncmp(option, "diffset=", 8) == 0 )
               {
                  if( diffsetname == 0 )
                  {
                     char* filename = &option[8];
                     diffsetname = strcpy(new char[strlen(filename) + 1], filename);
                  }
               }
               // --readmode=<value> : choose reading mode for <lpfile> (0* - floating-point, 1 - rational)
               else if( strncmp(option, "readmode=", 9) == 0 )
               {
                  if( !soplex.setIntParam(SoPlex::READMODE, option[9] - '0') )
                  {
                     freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
                     printUsage(argv, optidx);
                  }
                  // if the LP is parsed rationally and might be solved rationally, we choose automatic syncmode such that
                  // the rational LP is kept after reading
                  else if( soplex.intParam(SoPlex::READMODE) == SoPlex::READMODE_RATIONAL
                     && soplex.intParam(SoPlex::SOLVEMODE) != SoPlex::SOLVEMODE_REAL )
                  {
                     soplex.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
                  }
               }
               // --solvemode=<value> : choose solving mode (0* - floating-point solve, 1 - auto, 2 - force iterative refinement)
               else if( strncmp(option, "solvemode=", 10) == 0 )
               {
                  if( !soplex.setIntParam(SoPlex::SOLVEMODE, option[10] - '0') )
                  {
                     freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
                     printUsage(argv, optidx);
                  }
                  // if the LP is parsed rationally and might be solved rationally, we choose automatic syncmode such that
                  // the rational LP is kept after reading
                  else if( soplex.intParam(SoPlex::READMODE) == SoPlex::READMODE_RATIONAL
                     && soplex.intParam(SoPlex::SOLVEMODE) != SoPlex::SOLVEMODE_REAL )
                  {
                     soplex.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
                  }
               }
               // --<type>:<name>=<val> :  change parameter value using syntax of settings file entries
               else if( !soplex.parseSettingsString(option) )
               {
                  freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
                  printUsage(argv, optidx);
               }
               break;
            }

         case 't' :
            // -t<s> : set time limit to <s> seconds
            if( !soplex.setRealParam(SoPlex::TIMELIMIT, atoi(&option[2])) )
            {
               freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
               printUsage(argv, optidx);
            }
            break;

         case 'i' :
            // -i<n> : set iteration limit to <n>
            if( !soplex.setIntParam(SoPlex::ITERLIMIT, atoi(&option[2])) )
            {
               freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
               printUsage(argv, optidx);
            }
            break;

         case 'f' :
            // -f<eps> : set primal feasibility tolerance to <eps>
            if( !soplex.setRationalParam(SoPlex::FEASTOL, atof(&option[2])) )
            {
               freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
               printUsage(argv, optidx);
            }
            break;

         case 'o' :
            // -o<eps> : set dual feasibility (optimality) tolerance to <eps>
            if( !soplex.setRationalParam(SoPlex::OPTTOL, atof(&option[2])) )
            {
               freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
               printUsage(argv, optidx);
            }
            break;

         case 's' :
            // -s<value> : choose simplifier/presolver (0 - off, 1* - auto)
            if( !soplex.setIntParam(SoPlex::SIMPLIFIER, option[2] - '0') )
            {
               freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
               printUsage(argv, optidx);
            }
            break;

         case 'g' :
            // -g<value> : choose scaling (0 - off, 1 - uni-equilibrium, 2* - bi-equilibrium, 3 - geometric, 4 - iterated geometric)
            if( !soplex.setIntParam(SoPlex::SCALER, option[2] - '0') )
            {
               freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
               printUsage(argv, optidx);
            }
            break;

         case 'p' :
            // -p<value> : choose pricing (0* - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)
            if( !soplex.setIntParam(SoPlex::PRICER, option[2] - '0') )
            {
               freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
               printUsage(argv, optidx);
            }
            break;

         case 'v' :
            // -v<level> : set verbosity to <level> (0 - error, 3 - normal, 5 - high)
            if( !soplex.setIntParam(SoPlex::VERBOSITY, option[2] - '0') )
            {
               freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
               printUsage(argv, optidx);
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
               freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
               printUsage(argv, optidx);
            }
         }
      }

      // no LP file was given and no settings files are written
      if( optidx >= argc && savesetname == 0 && diffsetname == 0 )
      {
         freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
         printUsage(argv, 0);
      }

      // ensure that syncmode is not manual
      if( soplex.intParam(SoPlex::SYNCMODE) == SoPlex::SYNCMODE_MANUAL )
      {
         MSG_ERROR( spxout << "Error: manual synchronization is invalid on command line.  Change parameter int:syncmode.\n" );
         freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
         exit(1);
      }

      // save settings files
      if( savesetname != 0 )
      {
         MSG_INFO1( spxout << "Saving parameters to settings file <" << savesetname << "> . . .\n" );
         if( !soplex.saveSettingsFile(savesetname, false) )
         {
            MSG_ERROR( spxout << "Error writing parameters to file <" << savesetname << ">\n" );
         }
      }
      if( diffsetname != 0 )
      {
         MSG_INFO1( spxout << "Saving modified parameters to settings file <" << diffsetname << "> . . .\n" );
         if( !soplex.saveSettingsFile(diffsetname, true) )
         {
            MSG_ERROR( spxout << "Error writing modified parameters to file <" << diffsetname << ">\n" );
         }
      }

      // no LP file given: exit after saving settings
      if( optidx >= argc )
      {
         if( loadsetname != 0 || savesetname != 0 || diffsetname != 0 )
         {
            MSG_INFO1( spxout << "\n" );
         }

         freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
         exit(0);
      }

      // measure time for reading LP file and basis file
      readingTimer.start();

      // read LP from input file
      lpfilename = argv[optidx];
      MSG_INFO1( spxout << "Reading "
         << (soplex.intParam(SoPlex::READMODE) == SoPlex::READMODE_REAL ? "(real)" : "(rational)")
         << " LP file <" << lpfilename << "> . . .\n" );

      if( !soplex.readFile(lpfilename, &rownames, &colnames) )
      {
         MSG_ERROR( spxout << "Error while reading file <" << lpfilename << ">.\n" );
         freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
         exit(1);
      }

      // read basis file if specified
      if( readbasname != 0 )
      {
         MSG_INFO1( spxout << "Reading basis file <" << readbasname << "> . . . " );
         if( !soplex.readBasisFile(readbasname, &rownames, &colnames) )
         {
            MSG_ERROR( spxout << "Error while reading file <" << readbasname << ">.\n" );
            freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
            exit(1);
         }
      }

      readingTimer.stop();

      MSG_INFO1( std::streamsize prec = spxout.precision();
         spxout << "Reading took "
         << std::fixed << std::setprecision(2) << readingTimer.userTime()
         << std::scientific << std::setprecision(int(prec))
         << " seconds.\n\n" );

      MSG_INFO1( spxout << "LP has " << soplex.numRowsReal() << " rows "
         << soplex.numColsReal() << " columns and " << soplex.numNonzerosReal() << " nonzeros.\n\n" );

      // solve the LP
      soplex.solve();

      // print solution, check solution, and display statistics
      if( printPrimal )
      {
         ///@todo
      }

      if( printDual )
      {
         ///@todo
      }

      if( checkSol )
      {
         ///@todo
      }

      if( displayStatistics )
      {
         MSG_INFO1( spxout << "Statistics\n==========\n\n" );
         soplex.printStatistics(spxout.getStream(SPxOut::INFO1));
         MSG_INFO1( spxout << "\n" );
      }

      // write basis file if specified
      if( writebasname != 0 )
      {
         if( !soplex.hasBasis() )
         {
            MSG_WARNING( spxout << "No basis information available.  Could not write file <" << writebasname << ">\n\n" );
         }
         else if( !soplex.writeBasisFile(writebasname, &rownames, &colnames) )
         {
            MSG_ERROR( spxout << "Error while writing file <" << writebasname << ">.\n\n" );
            freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
            exit(1);
         }
         else
         {
            MSG_INFO1( spxout << "Written basis information to file <" << writebasname << ">.\n\n" );
         }
      }
   }
   catch( SPxException& x )
   {
      MSG_ERROR( spxout << "Exception caught: " << x.what() << "\n" );
      freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
      exit(1);
   }

   freeStrings(readbasname, writebasname, loadsetname, savesetname, diffsetname);
   exit(0);
}
