/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  args.hpp
 * @brief Code for argument parsing
 */

#ifndef _ARGS_HPP_
#define _ARGS_HPP_



/* Notes: If we need to add new parameters to things for Settings class. say,
   BoolParam:

          1. Add it to soplex.hpp under the name BoolParam default
          constructor.
          2. Add it to variable boolParam in this file.

   Similarly for other parameters.

   Unfortunately, for every new parameters in the settings class it has to be
   added at two places (soplex.cpp and soplex/args.hpp.) There is no way around
   this problem because, under current design the arguments are all parsed
   before a soplex object is created.

   For parameters that only affect the code, e.g., --checkfinal or -c (checks
   the final solution), this only has to be added inside the parseArgs function
   in this file. */

// Contains the list of all arguments
#include <initializer_list>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <initializer_list>
#include <functional>
#include <boost/program_options.hpp>
#include <boost/program_options/config.hpp> // For parse_config_file
#include <boost/program_options/errors.hpp>
#include <boost/exception/diagnostic_information.hpp>

#include <soplex/spxdefines.h>  // For access to some constants

#include <boost/multiprecision/number.hpp>

#ifdef SOPLEX_WITH_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

#ifdef SOPLEX_WITH_MPFR
// For multiple precision
#include <boost/multiprecision/mpfr.hpp>
#ifndef NDEBUG
#include "boost/multiprecision/debug_adaptor.hpp" // For debuging mpf numbers
#endif // NDEBUG
#endif // SOPLEX_WITH_MPFR
#ifdef SOPLEX_WITH_CPPMPF
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif  // SOPLEX_WITH_CPPMPF


namespace po = boost::program_options;

namespace soplex
{

// Runs SoPlex with parsed boost variables map, defined in soplexmain
template <class R>
int runSoPlex(const po::variables_map& vm);

namespace args
{

// A helper function to check if a command line argument lies in a range or in a
// list of values. Throws an exception if it doesn't. The functions return
// another function/lambda Will be used during the vm.notify()

// Returns a function that checks if the val lies in [min, max] TODO: If we
// have c++14/c++17, we can replace all the T with auto or use a template and
// put this inside the parseArgs. This also means that we can get rid of the
// ugly std::function part
template <typename T>
auto checkRange(const T& min, const T& max, const std::string& str) -> std::function<void(T)>
{
   // interesting remark: If you do &min, &max, &str instead of str below,
   // there will be a bug in some cases. Seems like the code is not
   // guaranteed to work according to the standard? Relevant SO:
   // https://stackoverflow.com/a/21443273/4223038
   return [min, max, str](const T & val)
   {
      if(val < min || val > max)
      {
         throw po::validation_error(po::validation_error::invalid_option_value,
                                    str + ", value=" + std::to_string(val));
      }
   };
}
} // namespace args ends here


// Parses the command line arguments
inline auto parseArgsAndRun(int argc, char* argv[]) -> int
{
   int solvemode = 1;
   unsigned int precision = 50;

   // a special case for working with ./soplex file.mps, i.e., without
   // explicitly doing ./soplex --lpfile=file.mps
   po::positional_options_description p;
   p.add("lpfile", -1);

   // Define all the options
   po::options_description generic("generic options");
   po::options_description algo("algorithmic settings (default in brackets)");
   po::options_description lt("limits and tolerances");
   po::options_description display("display options");
   po::options_description general("general options");
   po::options_description boolParam("bool Parameters");
   po::options_description intParam("Integer parameters");
   po::options_description realParam("Real parameters");
   // Rational param is currently empty. Can be used in the future
   po::options_description rationalParam("Rational parameters");


   generic.add_options()
   ("help,h", "help")
   ("allparam",
    "Displays complete set of parameters available (e.g., intParams, realParams boolParams))")
   ("version", "version");

   // Problem: Usually the descriptors should be just readbas. But how do I
   // print something like readbas=<basfile>?
   general.add_options()
   ("lpfile", po::value<std::string>(), "the lp file, specifying \"--lpfile\" is optional")
   ("readbas", po::value<std::string>(),  "read starting basis from file")
   ("writebas", po::value<std::string>(), "write terminal basis to file")
   ("writefile", po::value<std::string>(),
    "write LP to file in LP or MPS format depending on extension")
   ("writedual", po::value<std::string>(),
    "write the dual LP to a file in LP or MPS formal depending on extension")
   ("<type>:<name>=<val>", "change parameter value using syntax of settings file entries")
   // this is more of a hack. It won't be used. The only purpose of the above
   // argument is to print the thing.
   ("uint:random_seed", po::value<unsigned int>(), "set the random seed")
   // uint:random_seed used to be parsed inside parseSettingsString().
   ("loadset", po::value<std::string>(),
    "load parameters from settings file (overruled by command line parameters")
   ("saveset", po::value<std::string>(), "save parameters to settings file")
   ("diffset", po::value<std::string>(), "save modified parameters to settings file")
   ("extsol", po::value<std::string>(), "external solution for soplex to use for validation");

   lt.add_options()
   ("time,t", po::value<int>(), "set time limit to n seconds")
   ("iterlimit,i", po::value<int>(), "set iteration limit to n")
   ("primfeastol,f", po::value<Real>(), "set primal feasibility tolerance to double")
   ("dualfeastol,o", po::value<Real>(), "set dual feasibility (optimality) tolerance to")
   ("valtol,l", po::value<Real>(), "set validation tolerance to whatever");

   // Variables will contain value specified by user or the default

   algo.add_options()
   ("readmode", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 1, "readmode")),
    "choose reading mode for <lpfile> (0 - floating-point, 1 - rational)")
   ("solvemode", po::value<int>(&solvemode)->default_value(1)->notifier(args::checkRange(0, 4,
         "solvemode")),
    "choose solving mode (0 - floating-point, 1 - auto, 2 - force iterative refinement, 3 - multi-precision, 4 - Quadruple-precision (128-bit))")
   ("simplifier,s", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 3, "simplifier")),
    "choose simplifier/presolver (0 - off, 1 - auto)")
   ("scaler,g", po::value<int>()->default_value(2)->notifier(args::checkRange(0, 6, "scaler")),
    "choose scaling (0 - off, 1 - uni-equilibrium, 2 - bi-equilibrium, 3 - geometric, 4 - iterated geometric, 5 - least squares, 6 - geometric-equilibrium)")
   ("pricer,p", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 5, "pricer")),
    "choose pricing (0 - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)")
   ("ratiotester,r", po::value<int>()->default_value(3)->notifier(args::checkRange(0, 3,
         "ratiotester")),
    "choose ratio tester (0 - textbook, 1 - harris, 2 - fast, 3 - boundflipping)");

   display.add_options()
   ("verbosity,v", po::value<int>()->default_value(3)->notifier(args::checkRange(0, 5, "verbosity")),
    "set verbosity to <level> (0 - error, 3 - normal, 5 - high)")
   // Although the option says verbosity can be 0, 3, 5. In the program
   // verbosity 4 is also used.
   ("printprimal,x",  "print primal solution")
   ("printdualmult,y", "print dual multipliers")
   ("printprimratsol,X", "print primal solution in rational numbers")
   ("printdualmultrational,Y", "print dual multipliers in rational numbers")
   ("dispstat,q", "display detailed statistics")
   ("checkfinal,c", "perform final check of optimal solution in original problem");

   // These are parameters corresponding to the boolParam in
   // SoPlex<R>::Settings class. This doesn't have a notify function that
   // checks for the range because the boost can deal with it. It accepts the
   // parameters yes/no, on/off, 1/0 and true/false!
   boolParam.add_options()
   ("bool:lifting", po::value<bool>()->default_value(false),
    "should lifting be used to reduce range of nonzero matrix coefficients?")
   ("bool:eqtrans", po::value<bool>()->default_value(false),
    "should LP be transformed to equality form before a rational solve?")
   ("bool:testdualinf", po::value<bool>()->default_value(false),
    "should dual infeasibility be tested in order to try to return a dual solution even if primal infeasible?")
   ("bool:ratfac", po::value<bool>()->default_value(true),
    "should dual infeasibility be tested in order to try to return a dual solution even if primal infeasible?")
   ("bool:decompositiondualsimplex", po::value<bool>()->default_value(false),
    "should the decomposition based dual simplex be used to solve the LP?")
   ("bool:computedegen", po::value<bool>()->default_value(false),
    "should the degeneracy be computed for each basis?")
   ("bool:usecompdual", po::value<bool>()->default_value(false),
    "should the dual of the complementary problem be used in the decomposition simplex?")
   ("bool:explicitviol", po::value<bool>()->default_value(false),
    "Should violations of the original problem be explicitly computed in the decomposition simplex?")
   ("bool:acceptcycling", po::value<bool>()->default_value(false),
    "should cycling solutions be accepted during iterative refinement?")
   ("bool:ratrec", po::value<bool>()->default_value(true),
    "apply rational reconstruction after each iterative refinement?")
   ("bool:powerscaling", po::value<bool>()->default_value(true),
    "round scaling factors for iterative refinement to powers of two?")
   ("bool:ratfacjump", po::value<bool>()->default_value(false),
    "continue iterative refinement with exact basic solution if not optimal?")
   ("bool:rowboundflips", po::value<bool>()->default_value(false),
    "use bound flipping also for row representation?")
   ("bool:persistentscaling", po::value<bool>()->default_value(true),
    "should persistent scaling be used?")
   ("bool:fullperturbation", po::value<bool>()->default_value(false),
    "should perturbation be applied to the entire problem?")
   ("bool:ensureray", po::value<bool>()->default_value(false),
    "re-optimize the original problem to get a proof (ray) of infeasibility/unboundedness?")
   ("bool:forcebasic", po::value<bool>()->default_value(false),
    "try to enforce that the optimal solution is a basic solution");

   intParam.add_options()
   ("int:objsense", po::value<int>()->default_value(1)->notifier(args::checkRange(-1, 1,
         "int:objsense")),
    "objective sense (-1 - minimize, +1 - maximize)")
   ("int:representation", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 2,
         "int:representation")),
    "type of computational form (0 - auto, 1 - column representation, 2 - row representation)")
   ("int:algorithm", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 1,
         "int:algorithm")), "type of algorithm (0 - primal, 1 - dual)")
   ("int:factor_update_type", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 1,
         "int:factor_update_type")), "type of LU update (0 - eta update, 1 - Forrest-Tomlin update)")
   ("int:factor_update_max", po::value<int>()->default_value(0)->notifier(args::checkRange(0, INT_MAX,
         "int:factor_update_max")), "maximum number of LU updates without fresh factorization (0 - auto)")
   ("int:iterlimit", po::value<int>()->default_value(-1)->notifier(args::checkRange(-1, INT_MAX,
         "int:iterlimit")), "iteration limit (-1 - no limit)")
   ("int:reflimit", po::value<int>()->default_value(-1)->notifier(args::checkRange(-1, INT_MAX,
         "int:reflimit")), "refinement limit (-1 - no limit)")
   ("int:stallreflimit", po::value<int>()->default_value(-1)->notifier(args::checkRange(-1, INT_MAX,
         "int:stallreflimit")), "stalling refinement limit (-1 - no limit)")
   ("int:displayfreq", po::value<int>()->default_value(200)->notifier(args::checkRange(1, INT_MAX,
         "int:displayfreq")), "display frequency")
   ("int:verbosity", po::value<int>()->default_value(3)->notifier(args::checkRange(0, 5,
         "int:verbosity")),
    "verbosity level (0 - error, 1 - warning, 2 - debug, 3 - normal, 4 - high, 5 - full)")
   ("int:simplifier", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 1,
         "int:simplifier")), "simplifier (0 - off, 1 - auto)")
   ("int:scaler", po::value<int>()->default_value(2)->notifier(args::checkRange(0, 6, "int:scaler")),
    "scaling (0 - off, 1 - uni-equilibrium, 2 - bi-equilibrium, 3 - geometric, 4 - iterated geometric, 5 - least squares, 6 - geometric-equilibrium)")
   ("int:starter", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 3, "int:starter")),
    "crash basis generated when starting from scratch (0 - none, 1 - weight, 2 - sum, 3 - vector)")
   ("int:pricer", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 5, "int:pricer")),
    "pricing method (0 - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)")
   ("int:ratiotester", po::value<int>()->default_value(3)->notifier(args::checkRange(0, 3,
         "int:ratiotester")),
    "method for ratio test (0 - textbook, 1 - harris, 2 - fast, 3 - boundflipping)")
   ("int:syncmode", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 2,
         "int:syncmode")),
    "mode for synchronizing real and rational LP (0 - store only real LP, 1 - auto, 2 - manual)")
   ("int:readmode", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 1,
         "int:readmode")), "mode for reading LP files (0 - floating-point, 1 - rational)")
   ("int:solvemode", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 4,
         "int:solvemode")),
    "mode for iterative refinement strategy (0 - floating-point solve, 1 - auto, 2 - exact rational solve)")
   ("int:checkmode", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 2,
         "int:checkmode")),
    "mode for a posteriori feasibility checks (0 - floating-point check, 1 - auto, 2 - exact rational check)")
   ("int:timer", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 2, "int:timer")),
    "type of timer (1 - cputime, aka. usertime, 2 - wallclock time, 0 - no timing)")
   ("int:hyperpricing", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 2,
         "int:hyperpricing")), "mode for hyper sparse pricing (0 - off, 1 - auto, 2 - always)")
   ("int:ratfac_minstalls", po::value<int>()->default_value(2)->notifier(args::checkRange(0, INT_MAX,
         "int:ratfac_minstalls")),
    "minimum number of stalling refinements since last pivot to trigger rational factorization")
   ("int:leastsq_maxrounds", po::value<int>()->default_value(50)->notifier(args::checkRange(0, INT_MAX,
         "int:leastsq_maxrounds")),
    "maximum number of conjugate gradient iterations in least square scaling")
   ("int:solution_polishing", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 2,
         "int:solution_polishing")),
    "mode for solution polishing (0 - off, 1 - max basic slack, 2 - min basic slack)")
   ("int:decomp_iterlimit", po::value<int>()->default_value(100)->notifier(args::checkRange(1, INT_MAX,
         "int:decomp_iterlimit")),
    "the number of iterations before the decomposition simplex initialisation solve is terminated")
   ("int:decomp_maxaddedrows", po::value<int>()->default_value(500)->notifier(args::checkRange(1,
         INT_MAX, "int:decomp_maxaddedrows")),
    "maximum number of rows that are added to the reduced problem when using the decomposition based simplex")
   ("int:decomp_displayfreq", po::value<int>()->default_value(50)->notifier(args::checkRange(1,
         INT_MAX, "int:decomp_displayfreq")),
    "the frequency that the decomposition based simplex status output is displayed.")
   ("int:decomp_verbosity", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 5,
         "int:decomp_verbosity")),
    "the verbosity of decomposition based simplex (0 - error, 1 - warning, 2 - debug, 3 - normal, 4 - high, 5 - full).")
   ("int:printbasismetric", po::value<int>()->default_value(-1)->notifier(args::checkRange(-1, 3,
         "int:printbasismetric")),
    "print basis metric during the solve (-1 - off, 0 - condition estimate , 1 - trace, 2 - determinant, 3 - condition)")
   ("int:stattimer", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 2,
         "int:stattimer")),
    "measure for statistics, e.g. factorization time (0 - off, 1 - user time, 2 - wallclock time)");

   realParam.add_options()
   ("real:feastol", po::value<Real>()->default_value(1e-6)->notifier(args::checkRange(0.0, 0.1,
         "real:feastol")), "primal feasibility tolerance")
   ("real:opttol", po::value<Real>()->default_value(1e-6)->notifier(args::checkRange(0.0, 1.0,
         "real:opttol")), "dual feasibility tolerance")
   ("real:epsilon_zero", po::value<Real>()->default_value(DEFAULT_EPS_ZERO)->notifier(
       args::checkRange(0.0, 1.0, "real:epsilon_zero")), "general zero tolerance")
   ("real:epsilon_factorization",
    po::value<Real>()->default_value(DEFAULT_EPS_FACTOR)->notifier(args::checkRange(0.0, 1.0,
          "real:epsilon_factorization")), "zero tolerance used in factorization")
   ("real:epsilon_update", po::value<Real>()->default_value(DEFAULT_EPS_UPDATE)->notifier(
       args::checkRange(0.0, 1.0, "real:epsilon_update")),
    "zero tolerance used in update of the factorization")
   ("real:epsilon_pivot", po::value<Real>()->default_value(DEFAULT_EPS_PIVOT)->notifier(
       args::checkRange(0.0, 1.0, "real:epsilon_pivot")), "pivot zero tolerance used in factorization")
   ("real:infty", po::value<Real>()->default_value(DEFAULT_INFINITY)->notifier(args::checkRange(1e10,
         1e100, "real:infty")), "infinity threshold")
   ("real:timelimit", po::value<Real>()->default_value(DEFAULT_INFINITY)->notifier(args::checkRange(
            Real(0.0), DEFAULT_INFINITY, "real:timelimit")), "time limit in seconds")
   ("real:objlimit_lower", po::value<Real>()->default_value(-DEFAULT_INFINITY)->notifier(
       args::checkRange(-DEFAULT_INFINITY, DEFAULT_INFINITY, "real:objlimit_lower")),
    "lower limit on objective value")
   ("real:objlimit_upper", po::value<Real>()->default_value(DEFAULT_INFINITY)->notifier(
       args::checkRange(-DEFAULT_INFINITY, DEFAULT_INFINITY, "real:objlimit_upper")),
    "upper limit on objective value")
   ("real:fpfeastol", po::value<Real>()->default_value(1e-9)->notifier(args::checkRange(1e-12, 1.0,
         "real:fpfeastol")),
    "working tolerance for feasibility in floating-point solver during iterative refinement")
   ("real:fpopttol", po::value<Real>()->default_value(1e-9)->notifier(args::checkRange(1e-12, 1.0,
         "real:fpopttol")),
    "working tolerance for optimality in floating-point solver during iterative refinement")
   ("real:maxscaleincr", po::value<Real>()->default_value(1e25)->notifier(args::checkRange(Real(1.0),
         DEFAULT_INFINITY, "real:maxscaleincr")), "maximum increase of scaling factors between refinements")
   ("real:liftminval", po::value<Real>()->default_value(0.000976562)->notifier(args::checkRange(0.0,
         0.1, "real:liftminval")),
    "lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)")
   ("real:liftmaxval", po::value<Real>()->default_value(1024.0)->notifier(args::checkRange(Real(10.0),
         DEFAULT_INFINITY, "real:liftmaxval")),
    "lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)")
   ("real:sparsity_threshold", po::value<Real>()->default_value(0.6)->notifier(args::checkRange(0.0,
         1.0, "real:sparsity_threshold")),
    "sparse pricing threshold (#violations < dimension * SPARSITY_THRESHOLD activates sparse pricing)")
   ("real:representation_switch",
    po::value<Real>()->default_value(1.2)->notifier(args::checkRange(Real(0.0), DEFAULT_INFINITY,
          "real:representation_switch")),
    "threshold on number of rows vs. number of columns for switching from column to row representations in auto mode")
   ("real:ratrec_freq", po::value<Real>()->default_value(1.2)->notifier(args::checkRange(Real(1.0),
         DEFAULT_INFINITY, "real:ratrec_freq")),
    "geometric frequency at which to apply rational reconstruction")
   ("real:minred", po::value<Real>()->default_value(1e-4)->notifier(args::checkRange(0.0, 1.0,
         "real:minred")), "minimal reduction (sum of removed rows/cols) to continue simplification")
   ("real:refac_basis_nnz", po::value<Real>()->default_value(10.0)->notifier(args::checkRange(1.0,
         100.0, "real:refac_basis_nnz")),
    "refactor threshold for nonzeros in last factorized basis matrix compared to updated basis matrix")
   ("real:refac_update_fill", po::value<Real>()->default_value(5.0)->notifier(args::checkRange(1.0,
         100.0, "real:refac_update_fill")),
    "refactor threshold for fill-in in current factor update compared to fill-in in last factorization")
   ("real:refac_mem_factor", po::value<Real>()->default_value(1.5)->notifier(args::checkRange(1.0,
         10.0, "real:refac_mem_factor")),
    "refactor threshold for memory growth in factorization since last refactorization")
   ("real:leastsq_acrcy", po::value<Real>()->default_value(1000.0)->notifier(args::checkRange(Real(
            1.0), DEFAULT_INFINITY, "real:leastsq_acrcy")),
    "accuracy of conjugate gradient method in least squares scaling (higher value leads to more iterations)")
   ("real:obj_offset", po::value<Real>()->default_value(0.0)->notifier(args::checkRange(
            -DEFAULT_INFINITY, DEFAULT_INFINITY, "real:obj_offset")), "objective offset to be used")
   ("real:min_markowitz", po::value<Real>()->default_value(0.01)->notifier(args::checkRange(0.0001,
         0.9999, "real:min_markowitz")), "minimal Markowitz threshold in LU factorization");

#ifdef SOPLEX_WITH_MPFR
   po::options_description mpf("Multiprecision float solve");
   mpf.add_options()
   ("precision", po::value<unsigned int>(&precision)->default_value(50u),
    "Minimum precision (number of decimal digits) of mpf float. Only has effect if solvemode is set to 3");
#endif
#ifdef SOPLEX_WITH_CPPMPF
   po::options_description mpf("Multiprecision float solve");
   mpf.add_options()
   ("precision", po::value<unsigned int>(&precision)->default_value(50u),
    "Minimum precision (number of decimal digits) of cpp float. Only values 50, 100, 200 available. Compile with MPFR for other precisions. \n Only has effect if solvemode is set to 3");
#endif

   po::options_description allOpt("Allowed options");
   allOpt.add(generic).add(general).add(lt).add(algo).add(display).add(mpf).add(intParam).add(
      realParam).add(boolParam);

   // This will contain a subsection of the options, i.e., without intParam,
   // realParam, boolParam and rationalParam. Useful for printing a shorter help
   // message.
   po::options_description liteOpt("Allowed options");
   liteOpt.add(generic).add(general).add(lt).add(algo).add(display).add(mpf);

   try
   {
      // A hash map that takes a string to a boost::any type. The boost:any
      // type can be casted into the value we need.
      po::variables_map vm;
      po::store(po::command_line_parser(argc, argv).options(allOpt).positional(p).run(), vm);


      // If help is an option, print the help message and quit. This should
      // happen before notify
      if(vm.count("help"))
      {
         std::cout << liteOpt << "\n";
         return 0;
      }

      if(vm.count("allparam"))
      {
         std::cout << allOpt << "\n";
         return 0;
      }

      if(vm.count("version"))
      {
         SoPlexBase<Real> soplex;
         soplex.printVersion();
         return 0;
      }

      if(vm.count("loadset"))
      {
         const auto fname = vm["loadset"].as<std::string>();
         std::ifstream file(fname);

         if(file.fail())
         {
            BOOST_THROW_EXCEPTION(std::runtime_error("The settings file: " + fname + " does not exist"));
         }

         // According to the documentation if vm already has non-default values
         // it won't get replaced. This means the following:
         //
         // 1. soplex usually gets called by `./soplex --loadset=whatever.set
         // lpfile` So, one doesn't have to worry about losing the lpfile
         // information during the second call of store.
         //
         // 2. soplex is configured to run so that if you give additional
         // params and a settings file, like `./soplex --scaler=5
         // --loadset=whatever.set lpfile`, then the value of scaler must be 5
         // even if settings file has a different value for it.
         //
         // Both these problems won't occur.

         // Printing "Loading settings file..." happens inside the runSoPlex
         // function. A downside is that the reading time of the settings file
         // won't be part of the statistics. Hopefully this is okay since
         // settings files are small.
         po::store(po::parse_config_file(file, allOpt), vm);

      }

      // Relevant: https://stackoverflow.com/a/5517755/4223038 This will make
      // sure that all the required options are given, If not boost should
      // (throw) print an appropriate error. Also the notifier function while
      // making the argument makes sure that a condition is met, such as
      // whether the value given is in a certain range or inside an
      // initializer_list
      po::notify(vm);

      /* One of these has to be present for soplex to do something meaningful. */
      /* Help is displayed otherwise */
      if(!vm.count("lpfile") && !vm.count("diffset") && !vm.count("diffset") && !vm.count("loadset"))
      {
         std::cout << "Missing input file.\n\n";
         std::cout << liteOpt;
         return 0;
      }

      switch(solvemode)
      {
      case 0:                 // floating point
      case 1:                 // auto
      case 2:                 // iterative refinement
         runSoPlex<Real>(vm);
         break;

      case 3:                 // soplex mpf
         using namespace boost::multiprecision;

#if BOOST_VERSION < 107000
         std::cerr << "Error: Boost version too old." << std:: endl <<
                   "In order to use the multiprecision feature of SoPlex," <<
                   " Boost Version 1.70.0 or higher is required." << std::endl << \
                   "Included Boost version is " << BOOST_VERSION / 100000 << "."  // maj. version
                   << BOOST_VERSION / 100 % 1000 << "."  // min. version
                   << BOOST_VERSION % 100                // patch version;
                   << std::endl;
#else
#ifdef SOPLEX_WITH_MPFR

         // et_off means the expression templates options is turned off. TODO:
         // The documentation also mentions about static vs dynamic memory
         // allocation for the mpfr types. Is it relevant here? I probably also
         // need to have the mpfr_float_eto in the global soplex namespace
         using multiprecision = number<mpfr_float_backend<0>, et_off>;
         multiprecision::default_precision(precision);
         runSoPlex<multiprecision>(vm);
#endif  // SOPLEX_WITH_MPFR

#ifdef SOPLEX_WITH_CPPMPF
         // It seems that precision cannot be set on run time for cpp_float
         // backend for boost::number. So a precision of 50 decimal points is
         // set.
         using multiprecision1 = number<cpp_dec_float<50>, et_off>;
         using multiprecision2 = number<cpp_dec_float<100>, et_off>;
         using multiprecision3 = number<cpp_dec_float<200>, et_off>;

         if(precision <= 50)
            runSoPlex<multiprecision1>(vm);
         else if(precision <= 100)
            runSoPlex<multiprecision2>(vm);
         else
            runSoPlex<multiprecision3>(vm);

#endif  // SOPLEX_WITH_CPPMPF
#endif
         break;
#ifdef SOPLEX_WITH_FLOAT128

      case 4:                // quadprecision
#if BOOST_VERSION < 107000
         std::cerr << "Error: Boost version too old." << std:: endl <<
                   "In order to use the quadprecision feature of SoPlex," <<
                   " Boost Version 1.70.0 or higher is required." << std::endl << \
                   "Included Boost version is " << BOOST_VERSION / 100000 << "."  // maj. version
                   << BOOST_VERSION / 100 % 1000 << "."  // min. version
                   << BOOST_VERSION % 100                // patch version;
                   << std::endl;
#else
         using namespace boost::multiprecision;
         using Quad = boost::multiprecision::float128;
         runSoPlex<Quad>(vm);
#endif
         break;
#endif

      default:
         std::cerr << "Wrong value for the solve mode\n\n" << allOpt << "\n";
         return 0;
      }

   }
   catch(po::error& e)
   {
      // I think the whole verbosity thing won't apply here. So just a direct
      // call to std::cerr should be okay
      std::cerr << "Error in argument parsing: " << e.what() << "\n\n";
      // print the help message
      std::cout << liteOpt << "\n";
      return 1;
   }
   catch(...)
   {
      std::cerr << "Generic exception: " << boost::current_exception_diagnostic_information() << "\n";
      return 0;
   }

   return 0;
}

} // namespace soplex ends here

#endif // _ARGS_HPP_
