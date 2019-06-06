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

#include "boost/program_options.hpp"

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
   const char* usage =
      "general options:\n"
      "  --readbas=<basfile>    read starting basis from file\n"
      "  --writebas=<basfile>   write terminal basis to file\n"
      "  --writefile=<lpfile>   write LP to file in LP or MPS format depending on extension\n"
      "  --writedual=<lpfile>   write the dual LP to a file in LP or MPS formal depending on extension\n"
      "  --<type>:<name>=<val>  change parameter value using syntax of settings file entries\n"
      "  --loadset=<setfile>    load parameters from settings file (overruled by command line parameters)\n"
      "  --saveset=<setfile>    save parameters to settings file\n"
      "  --diffset=<setfile>    save modified parameters to settings file\n"
      "  --extsol=<value>       external solution for soplex to use for validation\n"
      "\n"
      "limits and tolerances:\n"
      "  -t<s>                  set time limit to <s> seconds\n"
      "  -i<n>                  set iteration limit to <n>\n"
      "  -f<eps>                set primal feasibility tolerance to <eps>\n"
      "  -o<eps>                set dual feasibility (optimality) tolerance to <eps>\n"
      "  -l<eps>                set validation tolerance to <eps>\n"
      "\n"
      "algorithmic settings (* indicates default):\n"
      "  --readmode=<value>     choose reading mode for <lpfile> (0* - floating-point, 1 - rational)\n"
      "  --solvemode=<value>    choose solving mode (0 - floating-point solve, 1* - auto, 2 - force iterative refinement)\n"
      "  -s<value>              choose simplifier/presolver (0 - off, 1* - auto)\n"
      "  -g<value>              choose scaling (0 - off, 1 - uni-equilibrium, 2* - bi-equilibrium, 3 - geometric, 4 - iterated geometric, 5 - least squares, 6 - geometric-equilibrium)\n"
      "  -p<value>              choose pricing (0* - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)\n"
      "  -r<value>              choose ratio tester (0 - textbook, 1 - harris, 2 - fast, 3* - boundflipping)\n"
      "\n"
      "display options:\n"
      "  -v<level>              set verbosity to <level> (0 - error, 3 - normal, 5 - high)\n"
      "  -x                     print primal solution\n"
      "  -y                     print dual multipliers\n"
      "  -X                     print primal solution in rational numbers\n"
      "  -Y                     print dual multipliers in rational numbers\n"
      "  -q                     display detailed statistics\n"
      "  -c                     perform final check of optimal solution in original problem\n"
      "\n";

   if(idx <= 0)
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
   if(s1 != 0)
   {
      delete [] s1;
      s1 = 0;
   }

   if(s2 != 0)
   {
      delete [] s2;
      s2 = 0;
   }

   if(s3 != 0)
   {
      delete [] s3;
      s3 = 0;
   }

   if(s4 != 0)
   {
      delete [] s4;
      s4 = 0;
   }

   if(s5 != 0)
   {
      delete [] s5;
      s5 = 0;
   }
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
void checkSolution(SoPlexBase<R>& soplex);

template <>
void checkSolution<Real>(SoPlexBase<Real>& soplex)
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
   using mpfr_float_100_eto = number<mpfr_float_backend<100>, et_off>;

   using mpfr_debug = number<debug_adaptor<mpfr_float_backend<100> >, et_off >;

   // The following won't compile. Because there is no conversion between mpq_t
   // Rational and cpp_float. Perhaps we need to change the Rational class

   // using cpp_float_50_eto = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<50>, boost::multiprecision::et_off>;

   // using cpp_float = boost::multiprecision::number<boost::multiprecision::cpp_

  //@todo need to implement the mpf part properly The arguments should be parsed
  // and the right template of runSoplex should be called

   // return runSoPlex<mpfr_debug>(argc, argv);
   // return runSoPlex<mpfr_float_100_eto>(argc, argv);
  return (runSoPlex<Real>(argc, argv)); // For the Real SoPlex

}

