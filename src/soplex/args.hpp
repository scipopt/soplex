// Contains the list of all arguments

namespace soplex
{

  auto parseArgs(int argc, char* argv[])
  {

    int verbosity;
    int readmode, solvemode, s, g, p, r;

    // Define all the options
    po::options_description generic("generic options");
    generic.add_options()
      ("help,h", "help")
      ("version", "version");

    // Problem: Usually the descriptors should be just readbas. But how do I
    // print something like readbas=<basfile>?
    po::options_description general("general options");
    general.add_options()
      ("readbas", po::value<std::vector<std::string> >(),  "read starting basis from file")
      ("writebas", po::value<std::vector<std::string> >(), "write terminal basis to file")
      ("writefile", po::value<std::vector<std::string> >(), "write LP to file in LP or MPS format depending on extension")
      ("writedual", po::value<std::vector<std::string> >(),  "write the dual LP to a file in LP or MPS formal depending on extension")
      ("<type>:<name>=<val>", "change parameter value using syntax of settings file entries") // TODO: How do I deal with this?
      ("loadset", po::value<std::vector<std::string> >(), "load parameters from settings file (overruled by command line parameters") // TODO: How do I deal with overruled?
      ("saveset", po::value<std::vector<std::string> >(), "save parameters to settings file")
      ("diffset", po::value<std::vector<std::string> >(), "save modified parameters to settings file")
      ("extsol", po::value<std::vector<std::string> >(), "external solution for soplex to use for validation");

    po::options_description lt("limits and tolerances");
    lt.add_options()
      ("time,t", po::value<int>(), "set time limit to n seconds")
      ("iterlimit,i", po::value<int>(), "set iteration limit to n")
      ("primfeastol,f", po::value<double>(), "set primal feasibility tolerance to double") // fix
      ("dualfeastol,o", po::value<double>(),"set dual feasibility (optimality) tolerance to") // vix
      ("valtol,l", po::value<double>(), "set validation tolerance to whatever"); // fix

    // Variables will contain value specified by user or the default

    po::options_description algo("algorithmic settings (default in brackets)");
    algo.add_options()
      ("readmode", po::value<int>(&readmode)->default_value(0), "choose reading mode for <lpfile> (0 - floating-point, 1 - rational)")
      ("solvemode", po::value<int>(&solvemode)->default_value(1), "choose solving mode (0 - floating-point solve, 1 - auto, 2 - force iterative
refinement, 3 - multi precision solve)")
      ("simplifier,s", po::value<int>(&s)->default_value(1), "choose simplifier/presolver (0 - off, 1 - auto)")
      ("scaler,g", po::value<int>(&g)->default_value(2),  "choose scaling (0 - off, 1 - uni-equilibrium, 2 - bi-equilibrium, 3 - geometric, 4 - iterated geometric, 5 - least squares, 6 - geometric-equilibrium)")
      ("pricer,p", po::value<int>(&p)->default_value(0), "choose pricing (0 - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)")
      ("ratiotester,r", po::value<int>(&r)->default_value(3), "choose ratio tester (0 - textbook, 1 - harris, 2 - fast, 3 - boundflipping)");

    po::options_description display("display options");
    display.add_options()
      ("verbosity,v", po::value<int>(&verbosity)->default_value(3), "set verbosity to <level> (0 - error, 3 - normal, 5 - high)") // fix the default
      ("printprimal,x",  "print primal solution")
      ("printdualmult,y", "print dual multipliers")
      ("printratsol,X", "print primal solution in rational numbers")
      ("printdualmultrational,Y", "print dual multipliers in rational numbers")
      ("detstat,q", "display detailed statistics")
      ("checkfinal,c", "perform final check of optimal solution in original problem");

    // These are parameters corresponding to the boolParam in SoPlex<R>::Settings
    // class. Ideally I want the parsing of arguments to finish before creating a
    // SoPlex Object. In the old SoPlex, half of the parsing happened outside
    // SoPlex and the rest inside the SoPlex.
    po::options_description boolParam("bool Parameters");
    param.add_options()
      ("bool:lifting", po::value<bool>()->default_value(false), "should lifting be used to reduce range of nonzero matrix coefficients?")
      ("bool:eqtrans", po::value<bool>()->default_value(false), "should LP be transformed to equality form before a rational solve?")
      ("testdualinf", po::value<bool>()->default_value(false), "should dual infeasibility be tested in order to try to return a dual solution even if primal infeasible?");

    po::options_description mpf("Multiprecision float solve");
    mpf.add_options()
      ("mpf", "Run templated multi-precision SoPlex") // This is redundant; there is the solvemode param
      ("precision", po::value<unsigned int>()->default_value(), "Minimum precision of mpf float");

    po::options_description allOpt("Allowed options");
    allOpt.add(generic).add(general).add(lt).add(algo).add(display);

    try
      {
        // This is like a hash table
        po::variables_map vm;
        // TODO: Maybe replace the parse_command_line with the other function?
        po::store(po::parse_command_line(argc, argv, allOpt), vm);

        po::notify(vm);

        // Idea. Parse all the program options without creating a SoPlex Object.
        // Check what is the solve mode and create a SoPlex Object accordingly
        // (i.e., Real V/S mpf)

        if(vm.count("help"))
          {
            std::cout<<allOpt<<"\n";
            return 0;
          }

        if(vm.count(",q"))         // This doesn't seem to be working. How do I deal with this?
          {
            std::cout<<"Display detailed statistics enabled\n";
          }

        if(vm.count("verbosity"))
          {
            std::cout<<"Verbosity is "<<v<<std::endl;
          }

        switch(solvemode)
          {
          case 0:                 // floating point
            // SoPlexBase<Real> soplex(); // Use the correct thing.
            std::cout<<"You are now going to do regular solve";
            break;
          case 1:                 // auto
            std::cout<<"You are now going to auto solve\n";
            break;
          case 2:                 // iterative refinement
            std::cout<<"You are now going to do iterative refinement\n";
            break;
          case 3:                 // soplex mpf
            // SoPlexBase<number<T> > soplex();
            std::cout<<"You are now running on mpf\n";
            break;
          default:
            std::cout<<"Wrong value for the solve mode\n\n"<<allOpt<<"\n";
            return 0;
          }

      }
    catch(std::exception& e)
      {
        std::cerr<<"error: "<<e.what()<<"\n\n";
        // print the help message
        std::cout<<allOpt<<"\n";
        return 1;
      }
    catch(...)
      {
        std::cerr<<"Exception of unknown type\n";
      }

    return 0;

  }
}
