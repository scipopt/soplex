// Contains the list of all arguments
#include <initializer_list>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <initializer_list>
#include <functional>
#include <boost/program_options.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/exception/diagnostic_information.hpp>

#include <soplex/spxdefines.h>  // For access to some constants

namespace po = boost::program_options;

namespace soplex
{

    // Runs SoPlex with parsed boost variables map, defined in soplexmain
  template <class R>
  int runSoPlex(const po::variables_map& vm);

  namespace args{
  // Checks if the val lies in [min, max] todo If we have c++14, we can replace
  // all the "int" with auto or use a template and put this inside the
  // parseArgs. This also means that we can get rid of the ugly std::function part
  template <typename T>
  auto checkRange (const T& min, const T& max, const std::string& str) -> std::function<void(T)>
                    {
                      return [&min, &max, &str](const T& val)
                             {
                               if(val < min || val > max)
                                 {
                                   std::cout<<min<<" "<<max<<"; "<<val<<" "<<str<<std::endl;
                                   throw po::validation_error(po::validation_error::invalid_option_value, str, std::to_string(val));
                                 }
                             };
                    };
  } // namespace args ends here


  // Parses the command line arguments
  auto parseArgs(int argc, char* argv[]) -> int
  {

    // Two helper functions to check if a command line argument lies in a range
    // or in a list of values. Throws an exception if it doesn't. The functions
    // return another function/lambda Will be used during the vm.notify()


    // Checks whether a value is inside a list and if not, it throws an error
    // todo If we have c++14, we can replace all the "int" with auto or use a
    // template, also std::cend() and std::cbegin()
    auto in = [](const std::initializer_list<int>& list, const std::string& str)
              {
                return [&list, &str](const int& val)
                       {
                         const auto lEnd = list.end();
                         const auto iter = std::find(list.begin(), list.end(), val);

                         if(iter == lEnd)            // meaning that val is not in the list
                           {
                             throw po::validation_error(po::validation_error::invalid_option_value, str, std::to_string(val));
                           }
                       };
              };


    int solvemode = 1;


    po::positional_options_description p;
    p.add("lpfile", -1);

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


    // Define all the options
    generic.add_options()
      ("help,h", "help")
      ("allparam", "Displays complete set of parameters available (e.g., intParams, realParams boolParams))")
      ("version", "version");

    // Problem: Usually the descriptors should be just readbas. But how do I
    // print something like readbas=<basfile>?
    general.add_options()
      ("lpfile", po::value<std::string>()->required(), "the lp file, specifying \"--lpfile\" is optional")
      ("readbas", po::value<std::string>(),  "read starting basis from file")
      ("writebas", po::value<std::string>(), "write terminal basis to file")
      ("writefile", po::value<std::string>(), "write LP to file in LP or MPS format depending on extension")
      ("writedual", po::value<std::string>(),  "write the dual LP to a file in LP or MPS formal depending on extension")
      ("<type>:<name>=<val>", "change parameter value using syntax of settings file entries") // TODO: How do I deal with this?
      ("loadset", po::value<std::string>(), "load parameters from settings file (overruled by command line parameters") // TODO: How do I deal with overruled?
      ("saveset", po::value<std::string>(), "save parameters to settings file")
      ("diffset", po::value<std::string>(), "save modified parameters to settings file")
      ("extsol", po::value<std::string>(), "external solution for soplex to use for validation");

    lt.add_options()
      ("time,t", po::value<int>(), "set time limit to n seconds")
      ("iterlimit,i", po::value<int>(), "set iteration limit to n")
      ("primfeastol,f", po::value<double>(), "set primal feasibility tolerance to double") // fix
      ("dualfeastol,o", po::value<double>(),"set dual feasibility (optimality) tolerance to") // vix
      ("valtol,l", po::value<double>(), "set validation tolerance to whatever"); // fix

    // Variables will contain value specified by user or the default

    algo.add_options()
      ("readmode", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 1, "readmode")),
       "choose reading mode for <lpfile> (0 - floating-point, 1 - rational)")
      ("solvemode", po::value<int>(&solvemode)->default_value(1)->notifier(args::checkRange(0, 3, "solvemode")),
       "choose solving mode (0 - floating-point solve, 1 - auto, 2 - force iterative refinement, 3 - multi precision solve)")
      ("simplifier,s", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 3, "simplifier")),
       "choose simplifier/presolver (0 - off, 1 - auto)")
      ("scaler,g", po::value<int>()->default_value(2)->notifier(args::checkRange(0, 6, "scaler")),
       "choose scaling (0 - off, 1 - uni-equilibrium, 2 - bi-equilibrium, 3 - geometric, 4 - iterated geometric, 5 - least squares, 6 - geometric-equilibrium)")
      ("pricer,p", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 5, "pricer")),
       "choose pricing (0 - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)")
      ("ratiotester,r", po::value<int>()->default_value(3)->notifier(args::checkRange(0, 3, "ratiotester")),
       "choose ratio tester (0 - textbook, 1 - harris, 2 - fast, 3 - boundflipping)");

    display.add_options()
      ("verbosity,v", po::value<int>()->default_value(3)->notifier(in({0, 3, 5}, "verbosity")), // TODO: Figure this out, needs to be replaced by a set or an initializer list
       "set verbosity to <level> (0 - error, 3 - normal, 5 - high)") // fix the default
      ("printprimal,x",  "print primal solution")
      ("printdualmult,y", "print dual multipliers")
      ("printratsol,X", "print primal solution in rational numbers")
      ("printdualmultrational,Y", "print dual multipliers in rational numbers")
      ("dispstat,q", "display detailed statistics")
      ("checkfinal,c", "perform final check of optimal solution in original problem");

    // These are parameters corresponding to the boolParam in SoPlex<R>::Settings
    // class. Ideally I want the parsing of arguments to finish before creating a
    // SoPlex Object. In the old SoPlex, half of the parsing happened outside
    // SoPlex and the rest inside the SoPlex.
    boolParam.add_options()
      ("bool:lifting", po::value<bool>()->default_value(false), "should lifting be used to reduce range of nonzero matrix coefficients?")
      ("bool:eqtrans", po::value<bool>()->default_value(false), "should LP be transformed to equality form before a rational solve?")
      ("bool:testdualinf", po::value<bool>()->default_value(false), "should dual infeasibility be tested in order to try to return a dual solution even if primal infeasible?")
      ("bool:ratfac", po::value<bool>()->default_value(false), "should dual infeasibility be tested in order to try to return a dual solution even if primal infeasible?")
      ("bool:decompositiondualsimplex", po::value<bool>()->default_value(false), "should the decomposition based dual simplex be used to solve the LP?")
      ("bool:computedegen", po::value<bool>()->default_value(false), "should the degeneracy be computed for each basis?")
      ("bool:usecompdual", po::value<bool>()->default_value(false), "should the dual of the complementary problem be used in the decomposition simplex?")
      ("bool:explicitviol", po::value<bool>()->default_value(false), "Should violations of the original problem be explicitly computed in the decomposition simplex?")
      ("bool:acceptcycling", po::value<bool>()->default_value(false), "should cycling solutions be accepted during iterative refinement?")
      ("bool:ratrec", po::value<bool>()->default_value(true), "apply rational reconstruction after each iterative refinement?")
      ("bool:powerscaling", po::value<bool>()->default_value(true), "round scaling factors for iterative refinement to powers of two?")
      ("bool:ratfacjump", po::value<bool>()->default_value(false), "continue iterative refinement with exact basic solution if not optimal?")
      ("bool:rowboundflips", po::value<bool>()->default_value(false), "use bound flipping also for row representation?")
      ("bool:persistentscaling", po::value<bool>()->default_value(true), "should persistent scaling be used?")
      ("bool:fullperturbation", po::value<bool>()->default_value(false), "should perturbation be applied to the entire problem?")
      ("bool:ensureray", po::value<bool>()->default_value(false), "re-optimize the original problem to get a proof (ray) of infeasibility/unboundedness?");

    intParam.add_options()
      ("int:objsense", po::value<int>()->default_value(1)->notifier(in({-1, 1}, "int:objsense")), "objective sense (-1 - minimize, +1 - maximize)")
      ("int:representation", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 2, "int:representation")), "type of computational form (0 - auto, 1 - column representation, 2 - row representation)")
      ("int:algorithm", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 1, "int:algorithm")), "type of algorithm (0 - primal, 1 - dual)")
      ("int:factor_update_type", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 1, "int:factor_update_type")), "type of LU update (0 - eta update, 1 - Forrest-Tomlin update)")
      ("int:factor_update_max", po::value<int>()->default_value(0)->notifier(args::checkRange(0, INT_MAX, "int:factor_update_max")), "maximum number of LU updates without fresh factorization (0 - auto)")
      ("int:iterlimit", po::value<int>()->default_value(-1)->notifier(args::checkRange(-1, INT_MAX, "int:iterlimit")), "iteration limit (-1 - no limit)")
      ("int:reflimit", po::value<int>()->default_value(-1)->notifier(args::checkRange(-1, INT_MAX, "int:reflimit")), "refinement limit (-1 - no limit)")
      ("int:stallreflimit", po::value<int>()->default_value(-1)->notifier(args::checkRange(-1, INT_MAX, "int:stallreflimit")), "stalling refinement limit (-1 - no limit)")
      ("int:displayfreq", po::value<int>()->default_value(200)->notifier(args::checkRange(1, INT_MAX, "int:displayfreq")), "display frequency")
      ("int:verbosity", po::value<int>()->default_value(3)->notifier(args::checkRange(0, 5, "int:verbosity")), "verbosity level (0 - error, 1 - warning, 2 - debug, 3 - normal, 4 - high, 5 - full)")
      ("int:simplifier", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 1, "int:simplifier")), "simplifier (0 - off, 1 - auto)")
      ("int:scaler", po::value<int>()->default_value(2)->notifier(args::checkRange(0, 6, "int:scaler")), "scaling (0 - off, 1 - uni-equilibrium, 2 - bi-equilibrium, 3 - geometric, 4 - iterated geometric, 5 - least squares, 6 - geometric-equilibrium)")
      ("int:starter", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 3, "int:starter")), "crash basis generated when starting from scratch (0 - none, 1 - weight, 2 - sum, 3 - vector)")
      ("int:pricer", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 5, "int:pricer")), "pricing method (0 - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)")
      ("int:ratiotester", po::value<int>()->default_value(3)->notifier(args::checkRange(0, 3, "int:ratiotester")), "method for ratio test (0 - textbook, 1 - harris, 2 - fast, 3 - boundflipping)")
      ("int:syncmode", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 2, "int:syncmode")), "mode for synchronizing real and rational LP (0 - store only real LP, 1 - auto, 2 - manual)")
      ("int:readmode", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 1, "int:readmode")), "mode for reading LP files (0 - floating-point, 1 - rational)")
      ("int:solvemode", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 2, "int:solvemode")), "mode for iterative refinement strategy (0 - floating-point solve, 1 - auto, 2 - exact rational solve)")
      ("int:checkmode", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 2, "int:checkmode")), "mode for a posteriori feasibility checks (0 - floating-point check, 1 - auto, 2 - exact rational check)")
      ("int:timer", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 2, "int:timer")), "type of timer (1 - cputime, aka. usertime, 2 - wallclock time, 0 - no timing)")
      ("int:hyperpricing", po::value<int>()->default_value(1)->notifier(args::checkRange(0, 2, "int:hyperpricing")), "mode for hyper sparse pricing (0 - off, 1 - auto, 2 - always)")
      ("int:ratfac_minstalls", po::value<int>()->default_value(2)->notifier(args::checkRange(0, INT_MAX, "int:ratfac_minstalls")), "minimum number of stalling refinements since last pivot to trigger rational factorization")
      ("int:leastsq_maxrounds", po::value<int>()->default_value(50)->notifier(args::checkRange(0, INT_MAX, "int:leastsq_maxrounds")), "maximum number of conjugate gradient iterations in least square scaling")
      ("int:solution_polishing", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 2, "int:solution_polishing")), "mode for solution polishing (0 - off, 1 - max basic slack, 2 - min basic slack)")
      ("int:decomp_iterlimit", po::value<int>()->default_value(100)->notifier(args::checkRange(1, INT_MAX, "int:decomp_iterlimit")), "the number of iterations before the decomposition simplex initialisation solve is terminated")
      ("int:decomp_maxaddedrows", po::value<int>()->default_value(500)->notifier(args::checkRange(1, INT_MAX, "int:decomp_maxaddedrows")), "maximum number of rows that are added to the reduced problem when using the decomposition based simplex")
      ("int:decomp_displayfreq", po::value<int>()->default_value(50)->notifier(args::checkRange(1, INT_MAX, "int:decomp_displayfreq")), "the frequency that the decomposition based simplex status output is displayed.")
      ("int:decomp_verbosity", po::value<int>()->default_value(0)->notifier(args::checkRange(0, 5, "int:decomp_verbosity")), "the verbosity of decomposition based simplex (0 - error, 1 - warning, 2 - debug, 3 - normal, 4 - high, 5 - full).")
      ("int:printbasismetric", po::value<int>()->default_value(-1)->notifier(args::checkRange(-1, 3, "int:printbasismetric")), "print basis metric during the solve (-1 - off, 0 - condition estimate , 1 - trace, 2 - determinant, 3 - condition)");

    realParam.add_options()
      ("real:feastol", po::value<double>()->default_value(1e-6)->notifier(args::checkRange(0.0, 0.1, "real:feastol")),"primal feasibility tolerance")
      ("real:opttol", po::value<double>()->default_value(1e-6)->notifier(args::checkRange(0.0, 1.0, "real:opttol")), "dual feasibility tolerance")
      ("real:epsilon_zero", po::value<double>()->default_value(DEFAULT_EPS_ZERO)->notifier(args::checkRange(0.0, 1.0, "real:epsilon_zero")), "general zero tolerance")
      ("real:epsilon_factorization", po::value<double>()->default_value(DEFAULT_EPS_FACTOR)->notifier(args::checkRange(0.0, 1.0, "real:epsilon_factorization")), "zero tolerance used in factorization")
      ("real:epsilon_update", po::value<double>()->default_value(DEFAULT_EPS_UPDATE)->notifier(args::checkRange(0.0, 1.0, "real:epsilon_update")), "zero tolerance used in update of the factorization")
      ("real:epsilon_pivot", po::value<double>()->default_value(DEFAULT_EPS_PIVOT)->notifier(args::checkRange(0.0, 1.0, "real:epsilon_pivot")), "pivot zero tolerance used in factorization")
      ("real:infty", po::value<double>()->default_value(DEFAULT_INFINITY)->notifier(args::checkRange(1e10, 1e100, "real:infty")), "infinity threshold")
      ("real:timelimit", po::value<double>()->default_value(DEFAULT_INFINITY)->notifier(args::checkRange(0.0, DEFAULT_INFINITY, "real:timelimit")), "time limit in seconds")
      ("real:objlimit_lower", po::value<double>()->default_value(-DEFAULT_INFINITY)->notifier(args::checkRange(-DEFAULT_INFINITY, DEFAULT_INFINITY, "real:objlimit_lower")), "lower limit on objective value")
      ("real:objlimit_upper", po::value<double>()->default_value(DEFAULT_INFINITY)->notifier(args::checkRange(-DEFAULT_INFINITY, DEFAULT_INFINITY, "real:objlimit_upper")), "upper limit on objective value")
      ("real:fpfeastol", po::value<double>()->default_value(1e-9)->notifier(args::checkRange(1e-12, 1.0, "real:fpfeastol")), "working tolerance for feasibility in floating-point solver during iterative refinement")
      ("real:fpopttol", po::value<double>()->default_value(1e-9)->notifier(args::checkRange(1e-12, 1.0, "real:fpopttol")), "working tolerance for optimality in floating-point solver during iterative refinement")
      ("real:maxscaleincr", po::value<double>()->default_value(1e25)->notifier(args::checkRange(1.0, DEFAULT_INFINITY, "real:maxscaleincr")), "maximum increase of scaling factors between refinements")
      ("real:liftminval", po::value<double>()->default_value(0.000976562)->notifier(args::checkRange(0.0, 0.1, "real:liftminval")), "lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)")
      ("real:liftmaxval", po::value<double>()->default_value(1024.0)->notifier(args::checkRange(10.0, DEFAULT_INFINITY, "real:liftmaxval")), "lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)")
      ("real:sparsity_threshold", po::value<double>()->default_value(0.6)->notifier(args::checkRange(0.0, 1.0, "real:sparsity_threshold")), "sparse pricing threshold (#violations < dimension * SPARSITY_THRESHOLD activates sparse pricing)")
      ("real:representation_switch", po::value<double>()->default_value(1.2)->notifier(args::checkRange(0.0, DEFAULT_INFINITY, "real:representation_switch")), "threshold on number of rows vs. number of columns for switching from column to row representations in auto mode")
      ("real:ratrec_freq", po::value<double>()->default_value(1.2)->notifier(args::checkRange(1.0, DEFAULT_INFINITY, "real:ratrec_freq")), "geometric frequency at which to apply rational reconstruction")
      ("real:minred", po::value<double>()->default_value(1e-4)->notifier(args::checkRange(0.0, 1.0, "real:minred")), "minimal reduction (sum of removed rows/cols) to continue simplification")
      ("real:refac_basis_nnz", po::value<double>()->default_value(10.0)->notifier(args::checkRange(1.0, 100.0, "real:refac_basis_nnz")), "refactor threshold for nonzeros in last factorized basis matrix compared to updated basis matrix")
  ("real:refac_update_fill", po::value<double>()->default_value(5.0)->notifier(args::checkRange(1.0, 100.0, "real:refac_update_fill")), "refactor threshold for fill-in in current factor update compared to fill-in in last factorization")
  ("real:refac_mem_factor", po::value<double>()->default_value(1.5)->notifier(args::checkRange(1.0, 10.0, "real:refac_mem_factor")), "refactor threshold for memory growth in factorization since last refactorization")
  ("real:leastsq_acrcy", po::value<double>()->default_value(1000.0)->notifier(args::checkRange(1.0, DEFAULT_INFINITY, "real:leastsq_acrcy")), "accuracy of conjugate gradient method in least squares scaling (higher value leads to more iterations)")
  ("real:obj_offset", po::value<double>()->default_value(0.0)->notifier(args::checkRange(-DEFAULT_INFINITY, DEFAULT_INFINITY, "real:obj_offset")), "objective offset to be used");

  po::options_description mpf("Multiprecision float solve");
  mpf.add_options()
  ("mpf", "Run templated multi-precision SoPlex") // This is redundant; there is the solvemode param
  ("precision", po::value<unsigned int>()->default_value(100), "Minimum precision of mpf float");

  po::options_description allOpt("Allowed options");
  allOpt.add(generic).add(general).add(lt).add(algo).add(display).add(intParam).add(realParam).add(boolParam);

  // This will contain a subsection of the options, i.e., without intParam, realParam, boolParam and rationalParam
  po::options_description liteOpt("Allowed options");
  liteOpt.add(generic).add(general).add(lt).add(algo).add(display);

  try
    {
      // A hash map that takes a string to a boost::any type. The boost:any
      // type can be casted into the value we need.
      po::variables_map vm;
      // TODO: Maybe replace the parse_command_line with the other function?
      po::store(po::command_line_parser(argc, argv).options(allOpt).positional(p).run(), vm);


      // If help is an option, print the help message and quit. This should
      // happen before notify
      if(vm.count("help"))
        {
          std::cout<<liteOpt<<"\n";
          return 0;
        }

      if(vm.count("allparam"))
        {
          std::cout<<allOpt<<"\n";
          return 0;
        }

      // Relevant: https://stackoverflow.com/a/5517755/4223038 This will make
      // sure that all the required options are given, If not boost should
      // (throw?) print an appropriate error. Also the notifier function while
      // making the argument makes sure that a condition is met, such as
      // whether the value given is in a certain range or inside an
      // initializer_list
      po::notify(vm);

      switch(solvemode)
        {
        case 0:                 // floating point
        case 1:                 // auto
        case 2:                 // iterative refinement
          runSoPlex<Real>(vm);
          break;
        case 3:                 // soplex mpf
          std::cout<<"You are now running on mpf\n"; // TODO: mpf
          break;
        default:
          std::cout<<"Wrong value for the solve mode\n\n"<<allOpt<<"\n";
          return 0;
        }

    }
  catch(po::error& e)
    {
      // TODO: How does SoPlex Handle the std::cerr?
      std::cerr<<"error: "<<e.what()<<"\n\n";
      // print the help message
      std::cout<<liteOpt<<"\n";
      return 1;
    }
  catch(...)
    {
      std::cerr<<"Unhandled exception\n"<<boost::current_exception_diagnostic_information()<<"\n";
      return 0;
    }

  return 0;
}

// TODO Look into this
// Checks whether a file exists
bool fileExists(const std::string& str) // Maybe throw an exception
{
  return static_cast<bool>(std::ifstream(str));
}


} // namespace soplex ends here
