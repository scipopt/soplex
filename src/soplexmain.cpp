/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2019 Konrad-Zuse-Zentrum                            */
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
#include "soplex/validation.h"
#include "soplex/statistics.h"
#include "soplex/args.hpp"      // For argument parsing

#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/debug_adaptor.hpp> // For debuging mpf numbers

#include <boost/program_options.hpp>
#include <exception>
#include <boost/exception/diagnostic_information.hpp>
#include <boost/exception/exception.hpp>

#ifdef SOPLEX_WITH_EGLIB
extern "C" {
#include "soplex/EGlib.h"
}
#else
#define EGlpNumStart() {}
#define EGlpNumClear() {}
#endif

using namespace soplex;

// function prototype
int main(int argc, char* argv[]);

// prints usage and command line options
static
void printUsage(const char* const argv[], int idx)
{

  // TODO: How is this supposed to work?
// #ifdef SOPLEX_WITH_ZLIB
//              << "  <lpfile>               linear program as .mps[.gz] or .lp[.gz] file\n\n"
// #else
//              << "  <lpfile>               linear program as .mps or .lp file\n\n"

// #endif
//              << usage;
}

/// performs external feasibility check with real type
///@todo implement external check; currently we use the internal methods for convenience
static
void checkSolutionReal(SoPlexBase<Real>& soplex)
{
   if(soplex.hasSol())
   {
      Real boundviol;
      Real rowviol;
      Real sumviol;

      if(soplex.getBoundViolation(boundviol, sumviol) && soplex.getRowViolation(rowviol, sumviol))
      {
         MSG_INFO1(soplex.spxout,
                   Real maxviol = boundviol > rowviol ? boundviol : rowviol;
                   bool feasible = (maxviol <= soplex.realParam(SoPlexBase<Real>::FEASTOL));
                   soplex.spxout << "Primal solution " << (feasible ? "feasible" : "infeasible")
                   << " in original problem (max. violation = " << std::scientific << maxviol
                   << std::setprecision(8) << std::fixed << ").\n"
                  );
      }
      else
      {
         MSG_INFO1(soplex.spxout, soplex.spxout << "Could not check primal solution.\n");
      }
   }
   else
   {
      MSG_INFO1(soplex.spxout, soplex.spxout << "No primal solution available.\n");
   }

   if(soplex.hasSol())
   {
      Real redcostviol;
      Real dualviol;
      Real sumviol;

      if(soplex.getRedCostViolation(redcostviol, sumviol) && soplex.getDualViolation(dualviol, sumviol))
      {
         MSG_INFO1(soplex.spxout,
                   Real maxviol = redcostviol > dualviol ? redcostviol : dualviol;
                   bool feasible = (maxviol <= soplex.realParam(SoPlexBase<Real>::OPTTOL));
                   soplex.spxout << "Dual solution " << (feasible ? "feasible" : "infeasible")
                   << " in original problem (max. violation = " << std::scientific << maxviol
                   << std::setprecision(8) << std::fixed << ").\n"
                  );
      }
      else
      {
         MSG_INFO1(soplex.spxout, soplex.spxout << "Could not check dual solution.\n");
      }
   }
   else
   {
      MSG_INFO1(soplex.spxout, soplex.spxout << "No dual solution available.\n");
   }
}

/// performs external feasibility check with rational type
///@todo implement external check; currently we use the internal methods for convenience
static
void checkSolutionRational(SoPlexBase<Real>& soplex)
{
   if(soplex.hasSol())
   {
      Rational boundviol;
      Rational rowviol;
      Rational sumviol;

      if(soplex.getBoundViolationRational(boundviol, sumviol)
            && soplex.getRowViolationRational(rowviol, sumviol))
      {
         MSG_INFO1(soplex.spxout,
                   Rational maxviol = boundviol > rowviol ? boundviol : rowviol;
                   bool feasible = (maxviol <= soplex.realParam(SoPlexBase<Real>::FEASTOL));
                   soplex.spxout << "Primal solution " << (feasible ? "feasible" : "infeasible") <<
                   " in original problem (max. violation = " << rationalToString(maxviol) << ").\n"
                  );
      }
      else
      {
         MSG_INFO1(soplex.spxout, soplex.spxout << "Could not check primal solution.\n");
      }
   }
   else
   {
      MSG_INFO1(soplex.spxout, soplex.spxout << "No primal solution available.\n");
   }

   if(soplex.hasSol())
   {
      Rational redcostviol;
      Rational dualviol;
      Rational sumviol;

      if(soplex.getRedCostViolationRational(redcostviol, sumviol)
            && soplex.getDualViolationRational(dualviol, sumviol))
      {
         MSG_INFO1(soplex.spxout,
                   Rational maxviol = redcostviol > dualviol ? redcostviol : dualviol;
                   bool feasible = (maxviol <= soplex.realParam(SoPlexBase<Real>::OPTTOL));
                   soplex.spxout << "Dual solution " << (feasible ? "feasible" : "infeasible") <<
                   " in original problem (max. violation = " << rationalToString(maxviol) << ").\n"
                  );
      }
      else
      {
         MSG_INFO1(soplex.spxout, soplex.spxout << "Could not check dual solution.\n");
      }
   }
   else
   {
      MSG_INFO1(soplex.spxout, soplex.spxout << "No dual solution available.\n");
   }
}

/// performs external feasibility check according to check mode
template <class R>

void checkSolution(SoPlexBase<R>& soplex)
{
   if(soplex.intParam(SoPlexBase<Real>::CHECKMODE) == SoPlexBase<Real>::CHECKMODE_RATIONAL
         || (soplex.intParam(SoPlexBase<Real>::CHECKMODE) == SoPlexBase<Real>::CHECKMODE_AUTO
             && soplex.intParam(SoPlexBase<Real>::READMODE) == SoPlexBase<Real>::READMODE_RATIONAL))
   {
      checkSolutionRational(soplex);
   }
   else
   {
      checkSolutionReal(soplex);
   }

   MSG_INFO1(soplex.spxout, soplex.spxout << "\n");
}

static
void printPrimalSolution(SoPlexBase<Real>& soplex, NameSet& colnames, NameSet& rownames,
                         bool real = true, bool rational = false)
{
   int printprec;
   int printwidth;
   printprec = (int) - log10(double(Param::epsilon()));
   printwidth = printprec + 10;

   if(real)
   {
      DVectorBase<Real> primal(soplex.numCols());

      if(soplex.getPrimalRay(primal))
      {
         MSG_INFO1(soplex.spxout, soplex.spxout << "\nPrimal ray (name, value):\n";)

         for(int i = 0; i < soplex.numCols(); ++i)
         {
            if(isNotZero(primal[i]))
            {
               MSG_INFO1(soplex.spxout, soplex.spxout << colnames[i] << "\t"
                         << std::setw(printwidth) << std::setprecision(printprec)
                         << primal[i] << std::endl;)
            }
         }

         MSG_INFO1(soplex.spxout, soplex.spxout << "All other entries are zero (within "
                   << std::setprecision(1) << std::scientific << Param::epsilon()
                   << std::setprecision(8) << std::fixed
                   << ")." << std::endl;)
      }
      else if(soplex.isPrimalFeasible() && soplex.getPrimal(primal))
      {
         int nNonzeros = 0;
         MSG_INFO1(soplex.spxout, soplex.spxout << "\nPrimal solution (name, value):\n";)

         for(int i = 0; i < soplex.numCols(); ++i)
         {
            if(isNotZero(primal[i]))
            {
               MSG_INFO1(soplex.spxout, soplex.spxout << colnames[i] << "\t"
                         << std::setw(printwidth) << std::setprecision(printprec)
                         << primal[i] << std::endl;)
               ++nNonzeros;
            }
         }

         MSG_INFO1(soplex.spxout, soplex.spxout << "All other variables are zero (within "
                   << std::setprecision(1) << std::scientific << Param::epsilon()
                   << std::setprecision(8) << std::fixed
                   << "). Solution has " << nNonzeros << " nonzero entries." << std::endl;)
      }
      else
         MSG_INFO1(soplex.spxout, soplex.spxout << "No primal information available.\n")
      }

   if(rational)
   {
      DVectorRational primal(soplex.numCols());

      if(soplex.getPrimalRayRational(primal))
      {
         MSG_INFO1(soplex.spxout, soplex.spxout << "\nPrimal ray (name, value):\n";)

         for(int i = 0; i < soplex.numCols(); ++i)
         {
            if(primal[i] != (Rational) 0)
            {
               MSG_INFO1(soplex.spxout, soplex.spxout << colnames[i] << "\t"
                         << std::setw(printwidth) << std::setprecision(printprec)
                         << primal[i] << std::endl;)
            }
         }

         MSG_INFO1(soplex.spxout, soplex.spxout << "All other entries are zero." << std::endl;)
      }

      if(soplex.isPrimalFeasible() && soplex.getPrimalRational(primal))
      {
         int nNonzeros = 0;
         MSG_INFO1(soplex.spxout, soplex.spxout << "\nPrimal solution (name, value):\n";)

         for(int i = 0; i < soplex.numColsRational(); ++i)
         {
            if(primal[i] != (Rational) 0)
            {
               MSG_INFO1(soplex.spxout, soplex.spxout << colnames[i] << "\t" << primal[i] << std::endl;)
               ++nNonzeros;
            }
         }

         MSG_INFO1(soplex.spxout, soplex.spxout << "All other variables are zero. Solution has "
                   << nNonzeros << " nonzero entries." << std::endl;)
      }
      else
         MSG_INFO1(soplex.spxout, soplex.spxout << "No primal (rational) solution available.\n")

      }
}

template <class R>
static
void printDualSolution(SoPlexBase<R>& soplex, NameSet& colnames, NameSet& rownames,
                       bool real = true, bool rational = false)
{
   int printprec;
   int printwidth;
   printprec = (int) - log10(double(Param::epsilon()));
   printwidth = printprec + 10;

   if(real)
   {
      DVector dual(soplex.numRows());

      if(soplex.getDualFarkas(dual))
      {
         MSG_INFO1(soplex.spxout, soplex.spxout << "\nDual ray (name, value):\n";)

         for(int i = 0; i < soplex.numRows(); ++i)
         {
            if(isNotZero(dual[i]))
            {
               MSG_INFO1(soplex.spxout, soplex.spxout << rownames[i] << "\t"
                         << std::setw(printwidth) << std::setprecision(printprec)
                         << dual[i] << std::endl;)
            }
         }

         MSG_INFO1(soplex.spxout, soplex.spxout << "All other entries are zero (within "
                   << std::setprecision(1) << std::scientific << Param::epsilon()
                   << std::setprecision(8) << std::fixed << ")." << std::endl;)
      }
      else if(soplex.isDualFeasible() && soplex.getDual(dual))
      {
         MSG_INFO1(soplex.spxout, soplex.spxout << "\nDual solution (name, value):\n";)

         for(int i = 0; i < soplex.numRows(); ++i)
         {
            if(isNotZero(dual[i]))
            {
               MSG_INFO1(soplex.spxout, soplex.spxout << rownames[i] << "\t"
                         << std::setw(printwidth) << std::setprecision(printprec)
                         << dual[i] << std::endl;)
            }
         }

         MSG_INFO1(soplex.spxout, soplex.spxout << "All other dual values are zero (within "
                   << std::setprecision(1) << std::scientific << Param::epsilon()
                   << std::setprecision(8) << std::fixed << ")." << std::endl;)

         DVector redcost(soplex.numCols());

         if(soplex.getRedCost(redcost))
         {
            MSG_INFO1(soplex.spxout, soplex.spxout << "\nReduced costs (name, value):\n";)

            for(int i = 0; i < soplex.numCols(); ++i)
            {
               if(isNotZero(redcost[i]))
               {
                  MSG_INFO1(soplex.spxout, soplex.spxout << colnames[i] << "\t"
                            << std::setw(printwidth) << std::setprecision(printprec)
                            << redcost[i] << std::endl;)
               }
            }

            MSG_INFO1(soplex.spxout, soplex.spxout << "All other reduced costs are zero (within "
                      << std::setprecision(1) << std::scientific << Param::epsilon()
                      << std::setprecision(8) << std::fixed << ")." << std::endl;)
         }
      }
      else
         MSG_INFO1(soplex.spxout, soplex.spxout << "No dual information available.\n")
      }

   if(rational)
   {
      DVectorRational dual(soplex.numRows());

      if(soplex.getDualFarkasRational(dual))
      {
         MSG_INFO1(soplex.spxout, soplex.spxout << "\nDual ray (name, value):\n";)

         for(int i = 0; i < soplex.numRows(); ++i)
         {
            if(dual[i] != (Rational) 0)
            {
               MSG_INFO1(soplex.spxout, soplex.spxout << rownames[i] << "\t"
                         << std::setw(printwidth)
                         << std::setprecision(printprec)
                         << dual[i] << std::endl;)
            }
         }

         MSG_INFO1(soplex.spxout, soplex.spxout << "All other entries are zero." << std::endl;)
      }

      if(soplex.isDualFeasible() && soplex.getDualRational(dual))
      {
         MSG_INFO1(soplex.spxout, soplex.spxout << "\nDual solution (name, value):\n";)

         for(int i = 0; i < soplex.numRowsRational(); ++i)
         {
            if(dual[i] != (Rational) 0)
               MSG_INFO1(soplex.spxout, soplex.spxout << rownames[i] << "\t" << dual[i] << std::endl;)
            }

         MSG_INFO1(soplex.spxout, soplex.spxout << "All other dual values are zero." << std::endl;)

         DVectorRational redcost(soplex.numCols());

         if(soplex.getRedCostRational(redcost))
         {
            MSG_INFO1(soplex.spxout, soplex.spxout << "\nReduced costs (name, value):\n";)

            for(int i = 0; i < soplex.numCols(); ++i)
            {
               if(redcost[i] != (Rational) 0)
                  MSG_INFO1(soplex.spxout, soplex.spxout << colnames[i] << "\t" << redcost[i] << std::endl;)
               }

            MSG_INFO1(soplex.spxout, soplex.spxout << "All other reduced costs are zero." << std::endl;)
         }
      }
      else
         MSG_INFO1(soplex.spxout, soplex.spxout << "No dual (rational) solution available.\n")
      }
}


/// runs SoPlexBase command line
int main(int argc, char* argv[])
{
   ///@todo the EGlib version info should be printed after the SoPlexBase version info
   // initialize EGlib's GMP memory management before any rational numbers are created
   EGlpNumStart();

   using namespace boost::multiprecision;

   // mpfr_flot_50 with expression template turned off

   // The following won't compile. Because there is no conversion between mpq_t
   // Rational and cpp_float. Perhaps we need to change the Rational class

   // using cpp_float_50_eto = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<50>, boost::multiprecision::et_off>;

   // using cpp_float = boost::multiprecision::number<boost::multiprecision::cpp_

  //@todo need to implement the mpf part properly The arguments should be parsed
  // and the right template of runSoplex should be called


   // return runSoPlex<mpfr_debug>(argc, argv);
   // return runSoPlex<mpfr_float_100_eto>(argc, argv);
  // return (runSoPlex<Real>(argc, argv)); // For the Real SoPlex

   return parseArgs(argc, argv);

}


namespace soplex{

// Runs SoPlex with the parsed boost variables map
  template <class R>
  int runSoPlex(const po::variables_map& vm)
  {
    SoPlexBase<R> soplex;
    Validation<R> validation;

   Timer* readingTime = nullptr;

   int optidx;

   // Stores different names
   std::string lpfilename, readbasname, writebasname, writefilename, writedualfilename, loadsetname, savesetname, diffsetname;

   bool printPrimal = false;
   bool printPrimalRational = false;
   bool printDual = false;
   bool printDualRational = false;
   bool displayStatistics = false;
   bool checkSol = false;

   int returnValue = 0;

   // For optional argument the following helper function would assign the correct value to the correct variable
   auto readIntoString = [&vm](std::string& var, const std::string str)
                         {
                           if(vm.count(str))
                             {
                               var = vm[str].as<std::string>();
                             }
                         };

   readIntoString(lpfilename, "lpfile");
   readIntoString(readbasname, "readbas");
   readIntoString(writebasname, "writebas");
   readIntoString(writefilename, "writefile");
   readIntoString(writedualfilename, "writedual");
   readIntoString(loadsetname, "loadset");
   readIntoString(savesetname, "saveset");
   readIntoString(diffsetname, "diffset");

   try
   {
      NameSet rownames;
      NameSet colnames;

      // create default timer (CPU time)
      readingTime = TimerFactory::createTimer(Timer::USER_TIME);

        soplex.printVersion();
        MSG_INFO1( soplex.spxout, soplex.spxout << SOPLEX_COPYRIGHT << std::endl << std::endl );

      // TODO: Read the settings file first. Remember that command line
      // arguments are meant to override the settings file. This would work,
      // because the stuff from settings will be done now? What about the
      // readbasis etc from before? TODO If settings provide readbasis then,
      // we'll have a problem.

      // We do a ranged for loop for every element in array intParam.
      // Afterwards, we do a search for it in the variables map, which is like a
      // std::map. I think the count part may be redundant because in the
      // args.hpp in the options_description::add_options(), I'm already
      // specifying the default value. But the soplex::settings object already
      // has this default value; one of them is redundant. The
      // variables_map::value() will return a boost::any object.

      for(int i = 0; i < SoPlexBase<R>::INTPARAM_COUNT; ++i)
        {
            const auto str = "int:" + soplex._currentSettings->intParam.name[i];
          if(vm.count(str))
            {
                soplex.parseSettingsString(str, vm[str].value());
            }
        }

      for(int i = 0; i < SoPlexBase<R>::REALPARAM_COUNT; ++i)
        {
            const auto str = "real:" + soplex._currentSettings->realParam.name[i];
          if(vm.count(str))
            {
                soplex.parseSettingsString(str, vm[str].value());
            }
        }

      for(int i = 0; i < SoPlexBase<R>::BOOLPARAM_COUNT; ++i)
        {
            const auto str = "int:" + soplex._currentSettings->boolParam.name[i];
          if(vm.count(str))
            {
                soplex.parseSettingsString(str, vm[str].value());
            }
        }
      // I didn't write the loop for rationalParams, because currently there are
      // no rationalParams.


      if(vm.count("readmode"))
        {
            soplex.setIntParam(soplex.READMODE, vm["readmode"].as<int>());
        }

      // --solvemode=<value> : choose solving mode (0* - floating-point solve, 1 - auto, 2 - force iterative refinement)
      if(vm.count("solvemode"))
        {
            soplex.setIntParam(soplex.SOLVEMODE, vm["solvemode"].as<int>());

          // if the LP is parsed rationally and might be solved rationally, we choose automatic syncmode such that
          // the rational LP is kept after reading

          // TODO: Look at this if statement. It used to be an else if

            if( soplex.intParam(soplex.READMODE) == soplex.READMODE_RATIONAL
                && soplex.intParam(soplex.SOLVEMODE) != soplex.SOLVEMODE_REAL )
            {
                soplex.setIntParam(soplex.SYNCMODE, soplex.SYNCMODE_AUTO);
            }

        }

      // --extsol=<value> : external solution for soplex to use for validation
      if(vm.count("extsol"))
        {
          auto input = vm["extsol"].as<std::string>();

            validation.updateExternalSolution(input);
        }

      // settings file format arguments are handled at the beginning.

      // -t<s> : set time limit to <s> seconds
      if(vm.count("time"))
        {
            soplex.setRealParam(soplex.TIMELIMIT, vm["time"].as<int>());
        }

      // -i<n> : set iteration limit to <n>
      if(vm.count("iterlimit"))
        {
            soplex.setIntParam(soplex.ITERLIMIT, vm["iterlimit"].as<int>());
        }

      // -f<eps> : set primal feasibility tolerance to <eps>
      if(vm.count("primfeastol"))
        {
            soplex.setRealParam(soplex.FEASTOL, vm["primfeastol"].as<double>());
        }

      // -o<eps> : set dual feasibility (optimality) tolerance to <eps>
      if(vm.count("dualfeastol"))
        {
            soplex.setRealParam(soplex.OPTTOL, vm["dualfeastol"].as<double>());
        }

      // l<eps> : set validation tolerance to <eps>
      if(vm.count("valtol"))
        {
          auto str = vm["valtol"].as<std::string>();
            validation.updateValidationTolerance(str);
        }

      // -s<value> : choose simplifier/presolver (0 - off, 1* - auto)
      if(vm.count("simplifier"))
        {
            soplex.setIntParam(soplex.SIMPLIFIER, vm["simplifier"].as<int>());
        }

      // -g<value> : choose scaling (0 - off, 1 - uni-equilibrium, 2* - bi-equilibrium, 3 - geometric, 4 - iterated geometric,  5 - least squares, 6 - geometric-equilibrium)
      if(vm.count("scaler"))
        {
            soplex.setIntParam(soplex.SCALER, vm["scaler"].as<int>());
        }

      // -p<value> : choose pricing (0* - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)
      if(vm.count("pricer"))
     {
            soplex.setIntParam(soplex.PRICER, vm["pricer"].as<int>());
     }

      // -r<value> : choose ratio tester (0 - textbook, 1 - harris, 2* - fast, 3 - boundflipping)
      if(vm.count("ratiotester"))
        {
            soplex.setIntParam(soplex.RATIOTESTER, vm["ratiotester"].as<int>());
        }

      // -v<level> : set verbosity to <level> (0 - error, 3 - normal, 5 - high)
      if(vm.count("verbosity"))
        {
            soplex.setIntParam(soplex.VERBOSITY, vm["verbosity"].as<int>());
        }

      // -x : print primal solution
      if(vm.count("printprimal"))
        {
          printPrimal = true;
        }

      // -X : print primal solution with rationals
      if(vm.count("printprimratsol"))
        {
          printPrimalRational = true;
        }

      // -y : print dual multipliers
      if(vm.count("printdualmult"))
        {
          printDual = true;
        }

      // -Y : print dual multipliers with rationals
      if(vm.count("printdualmultrational"))
        {
          printDualRational = true;
        }

      // -q : display detailed statistics
      if(vm.count("dispstat"))
        {
          displayStatistics = true;
        }

      // -c : perform final check of optimal solution in original problem
      if(vm.count("checkfinal"))
        {
          checkSol = true;
        }

      if(vm.count("loadset"))
        {
          // This is more of a hack. The actual reading happens inside args.hpp
          // This is here because 1. the spxout object is not available there.
          // 2. It has to be displayed after the soplex intro banner
            MSG_INFO1(soplex.spxout, soplex.spxout << "Loading settings file <" << vm["loadset"].as<std::string>() << "> . . .\n");
        }

        // For random seed
        if(vm.count("uint:random_seed"))
          {
            soplex.setRandomSeed(vm["uint:random_seed"].as<unsigned int>());
          }

        MSG_INFO1( soplex.spxout, soplex.printUserSettings(); );

      // TODO: How is the following code supposed to work?
      // no LP file was given and no settings files are written
      if( lpfilename.empty() && savesetname.empty() && diffsetname.empty())
      {
        BOOST_THROW_EXCEPTION(std::invalid_argument("Empty lpfilename, save setings name and diff settings name"));
      }

      // ensure that syncmode is not manual
        if( soplex.intParam(soplex.SYNCMODE) == soplex.SYNCMODE_MANUAL )
      {
        BOOST_THROW_EXCEPTION(std::invalid_argument("Error: manual synchronization is invalid on command line.  Change parameter int:syncmode"));
      }

      // save settings files
      if(!savesetname.empty())
      {
            MSG_INFO1( soplex.spxout, soplex.spxout << "Saving parameters to settings file <" << savesetname << "> . . .\n" );
            if( !soplex.saveSettingsFile(savesetname.c_str(), false) )
         {
            MSG_ERROR( std::cerr << "Error writing parameters to file <" << savesetname << ">\n" );
         }
      }
      if(!diffsetname.empty())
      {
            MSG_INFO1( soplex.spxout, soplex.spxout << "Saving modified parameters to settings file <" << diffsetname << "> . . .\n" );
            if( !soplex.saveSettingsFile(diffsetname.c_str(), true) )
         {
            MSG_ERROR( std::cerr << "Error writing modified parameters to file <" << diffsetname << ">\n" );
         }
      }

      // no LP file given: exit after saving settings
      if(lpfilename.empty())
      {
        if(!loadsetname.empty() || !savesetname.empty() || !diffsetname.empty())
         {
                MSG_INFO1( soplex.spxout, soplex.spxout << "\n" );
         }

        BOOST_THROW_EXCEPTION(std::invalid_argument("No lipfile, settings file, save settings file or diff settings file"));
      }

      // measure time for reading LP file and basis file
      readingTime->start();

      // if the LP is parsed rationally and might be solved rationally, we choose automatic syncmode such that
      // the rational LP is kept after reading
        if( soplex.intParam(soplex.READMODE) == soplex.READMODE_RATIONAL
            && soplex.intParam(soplex.SOLVEMODE) != soplex.SOLVEMODE_REAL )
      {
            soplex.setIntParam(soplex.SYNCMODE, soplex.SYNCMODE_AUTO);
      }

      // read LP from input file
        MSG_INFO1( soplex.spxout, soplex.spxout << "Reading "
                   << (soplex.intParam(soplex.READMODE) == soplex.READMODE_REAL ? "(real)" : "(rational)")
         << " LP file <" << lpfilename << "> . . .\n" );

        if( !soplex.readFile(lpfilename.c_str(), &rownames, &colnames) )
      {
         BOOST_THROW_EXCEPTION(std::logic_error("Error while reading lpfile: " + lpfilename ));
      }

      // write LP if specified
      if(!writefilename.empty())
      {
            if( !soplex.writeFile(writefilename.c_str(), &rownames, &colnames) )
         {
            BOOST_THROW_EXCEPTION(std::runtime_error("Error in writing to file: " + writefilename));
         }
         else
         {
                MSG_INFO1( soplex.spxout, soplex.spxout << "Written LP to file <" << writefilename << ">.\n\n" );
         }
      }

      // write dual LP if specified
      if(!writedualfilename.empty())
      {
            if( !soplex.writeDualFileReal(writedualfilename.c_str(), &rownames, &colnames) )
         {
            BOOST_THROW_EXCEPTION(std::runtime_error("Error while writing dual file: " + writedualfilename));
         }
         else
         {
                MSG_INFO1( soplex.spxout, soplex.spxout << "Written dual LP to file <" << writedualfilename << ">.\n\n" );
         }
      }

      // read basis file if specified
      if(!readbasname.empty())
      {
            MSG_INFO1( soplex.spxout, soplex.spxout << "Reading basis file <" << readbasname << "> . . . " );
            if( !soplex.readBasisFile(readbasname.c_str(), &rownames, &colnames) )
         {
           BOOST_THROW_EXCEPTION(std::runtime_error("Error while reading file: " + readbasname));
         }
      }

      readingTime->stop();

        MSG_INFO1( soplex.spxout,
                   std::streamsize prec = soplex.spxout.precision();
                   soplex.spxout << "Reading took "
         << std::fixed << std::setprecision(2) << readingTime->time()
         << std::scientific << std::setprecision(int(prec))
         << " seconds.\n\n" );

        MSG_INFO1( soplex.spxout, soplex.spxout << "LP has " << soplex.numRows() << " rows "
                   << soplex.numCols() << " columns and " << soplex.numNonzeros() << " nonzeros.\n\n" );

      // solve the LP
        soplex.optimize();

      // print solution, check solution, and display statistics
        printPrimalSolution(soplex, colnames, rownames, printPrimal, printPrimalRational);
        printDualSolution(soplex, colnames, rownames, printDual, printDualRational);

      if( checkSol )
          checkSolution<R>(soplex); // The type needs to get fixed here

      if( displayStatistics )
      {
            MSG_INFO1( soplex.spxout, soplex.spxout << "Statistics\n==========\n\n" );
            soplex.printStatistics(soplex.spxout.getStream(SPxOut::INFO1));
      }

        if(validation.validate)
          validation.validateSolveReal(soplex);

      // write basis file if specified
      if(!writebasname.empty())
      {
            if( !soplex.hasBasis() )
         {
                MSG_WARNING( soplex.spxout, soplex.spxout << "No basis information available.  Could not write file <" << writebasname << ">\n\n" );
         }
            else if( !soplex.writeBasisFile(writebasname.c_str(), &rownames, &colnames) )
         {
            BOOST_THROW_EXCEPTION(std::runtime_error("Error while writing file: " + writebasname));
         }
         else
         {
                MSG_INFO1( soplex.spxout, soplex.spxout << "Written basis information to file <" << writebasname << ">.\n\n" );
         }
      }
   }
   catch( const SPxException& x ) // There could be an exception from boost?
   {
      MSG_ERROR( std::cerr << "Exception caught: " << x.what() << "\n" );
      returnValue = 1;
   }
   catch(...)
     {
        std::cout<<"Generic exception: "<<boost::current_exception_diagnostic_information()<<"\n";
       returnValue = 1;
     }


   // because EGlpNumClear() calls mpq_clear() for all mpq_t variables, we need to destroy all objects of class Rational
   // beforehand; hence all Rational objects and all data that uses Rational objects must be allocated dynamically via
   // spx_alloc() and freed here; disabling the list memory is crucial
 Rational::disableListMem();
 EGlpNumClear();
 if( nullptr != readingTime )
 {
    readingTime->~Timer();
    spx_free(readingTime);
 }

   return returnValue;
  }

} // namespace soplex ends here
// TODO: Is this required?

