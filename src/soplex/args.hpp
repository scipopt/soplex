// Contains the list of all arguments

namespace soplex
{
  namespace po = boost::program_options;

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
        // A hash map that takes a string to a boost::any type. The boost:any
        // type can be casted into the value we need.
        po::variables_map vm;
        // TODO: Maybe replace the parse_command_line with the other function?
        po::store(po::parse_command_line(argc, argv, allOpt), vm);


        // If help is an option, print the help message and quit
        if(vm.count("help"))
          {
            std::cout<<allOpt<<"\n";
            return 0;
          }

        // Relevant: https://stackoverflow.com/a/5517755/4223038
        // This will make sure that all the required options are given, If not
        // boost should (throw?) print an appropriate error
        po::notify(vm);

        switch(solvemode)
          {
          case 0:                 // floating point
          case 1:                 // auto
          case 2:                 // iterative refinement
            runSoPlex<Real>(vm);
            break;
          case 3:                 // soplex mpf
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


  // TODO Maybe I can have this function inside the soplexmain.cpp. And the args function would return a variables_map. I call the runSoPlex function with the variables map.
template <class R>
int runSoPlex(const po::variables_map& vm)
{
  SoPlexBase<R>* soplex = nullptr;

   Timer* readingTime = nullptr;
  Validation<R>* validation = nullptr;
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
                               // TODO Did I use the spxSnprintf correctly?
                               // Maybe the + 1 is wrong. Check out the syntax later
                               spxSnprintf(var, var.size() + 1, "%s", filename);
                             }
                         };

   // TODO: Figure out how to deal with lpfilename; it was mentioned in the documentation


   readIntoString(readbasname, "readbas");
   readIntoString(writebasname, "writebas");
   readIntoString(writefilename, "writefile");
   readIntoString(writedualfilename, "writedual");
   readIntoString(loadsetname, "loadset");
   readIntoString(savesetname, "saveset");
   readIntoString(diffsetname, "diffset");

   // TODO: What about extsol?


   try
   {
      NameSet rownames;
      NameSet colnames;

      // create default timer (CPU time)
      readingTime = TimerFactory::createTimer(Timer::USER_TIME);
      soplex = nullptr;
      spx_alloc(soplex);
      new (soplex) SoPlexBase<R>();

      soplex->printVersion();
      MSG_INFO1( soplex->spxout, soplex->spxout << SOPLEX_COPYRIGHT << std::endl << std::endl );

      validation = nullptr;
      spx_alloc(validation);
      new (validation) Validation<R>();


      // TODO: Read the settings file first. Remember that command line
      // arguments are meant to override the settings file. This would work,
      // because the stuff from settings will be done now? What about the
      // readbasis etc from before? TODO If settings provide readbasis then,
      // we'll have a problem.

      if(vm.count("readmode"))
         {
          // TODO how do we deal with the return value? I think the function
          // should thrown an exception if the value doesn't fit in the required
          // range.
          // TODO https://stackoverflow.com/q/25548090/4223038
          // The above stackoverflow link talks about the problem
          soplex->setIntParam(soplex->READMODE, vm["readmode"].as<int>());
         }

               // --solvemode=<value> : choose solving mode (0* - floating-point solve, 1 - auto, 2 - force iterative refinement)
      if(vm.count("solvemode"))
                  {
          soplex->setIntParam(soplex->SOLVEMODE, vm["solvemode"].as<int>());

                  // if the LP is parsed rationally and might be solved rationally, we choose automatic syncmode such that
                  // the rational LP is kept after reading
                  else if( soplex->intParam(soplex->READMODE) == soplex->READMODE_RATIONAL
                     && soplex->intParam(soplex->SOLVEMODE) != soplex->SOLVEMODE_REAL )
                  {
                     soplex->setIntParam(soplex->SYNCMODE, soplex->SYNCMODE_AUTO);
                  }

               }

               // --extsol=<value> : external solution for soplex to use for validation
      if(vm.count("extsol"))
                  {
          auto input = vm["extsol"].as<std::string>();

          validation.updateExternalSolution(input.c_str());
            }

      // TODO: Code for --type:name=<val>: How am I supposed to handle this?

            // -t<s> : set time limit to <s> seconds
      if(vm.count("time"))
            {
          soplex->setRealParam(soplex->TIMELIMIT, vm["time"].as<int>());
            }

            // -i<n> : set iteration limit to <n>
      if(vm.count("iterlimit"))
            {
          soplex->setRealParam(soplex->ITERLIMIT, vm["iterlimit"].as<int>());
            }

            // -f<eps> : set primal feasibility tolerance to <eps>
      if(vm.count("primfeastol"))
            {
          soplex->setRealParam(soplex->FEASTOL, vm["primfeastol"].as<double>());
            }

            // -o<eps> : set dual feasibility (optimality) tolerance to <eps>
      if(vm.count("dualfeastol"))
            {
          soplex->setRealParam(soplex->OPTTOL, vm["dualfeastol"].as<double>());
            }

            // l<eps> : set validation tolerance to <eps>
      if(vm.count("valtol"))
            {
          validation->updateValidationTolerance();
            }

            // -s<value> : choose simplifier/presolver (0 - off, 1* - auto)
      if(vm.count("simplifier"))
            {
          soplex->setIntParam(soplex->SIMPLIFIER, vm["simplifier"].as<int>());
            }

            // -g<value> : choose scaling (0 - off, 1 - uni-equilibrium, 2* - bi-equilibrium, 3 - geometric, 4 - iterated geometric,  5 - least squares, 6 - geometric-equilibrium)
      if(vm.count("scaler"))
            {
          soplex->setIntParam(soplex->PRICER, vm["scaler"]);
            }

            // -r<value> : choose ratio tester (0 - textbook, 1 - harris, 2* - fast, 3 - boundflipping)
      if(vm.count("ratiotester"))
            {
          soplex->setIntParam(soplex->RATIOTESTER, vm["rationtester"].as<int>());
            }

            // -v<level> : set verbosity to <level> (0 - error, 3 - normal, 5 - high)
      if(vm.count("verbosity"))
        {
          soplex->setIntParam(soplex->VERBOSITY, vm["verbosity"].as<int>());
        }

      // For the boolean variables printPrimal, printPrimalRational etc, they
      // will automatically get the values when the vm.notify() function gets
      // called

      // TODO what's the deal with the following code from the above stuff? case If --help or -h is called, then

      // 'h' : if( !soplex->saveSettingsFile(0, false) ) { MSG_ERROR( std::cerr
      // << "Error printing parameters\n" ); }


      // // read arguments from command line
      // for( optidx = 1; optidx < argc; optidx++ )
      // {
      //    char* option = argv[optidx];

      //    // we reached <lpfile>
      //    if( option[0] != '-' )
      //    {
      //       lpfilename = argv[optidx];
      //       continue;
      //    }

      //    // option string must start with '-', must contain at least two characters, and exactly two characters if and
      //    // only if it is -x, -y, -q, or -c
      //    if( option[0] != '-' || option[1] == '\0'
      //       || ((option[2] == '\0') != (option[1] == 'x' || option[1] == 'X' || option[1] == 'y' || option[1] == 'Y' || option[1] == 'q' || option[1] == 'c')) )
      //    {
      //       printUsage(argv, optidx);
      //       returnValue = 1;
      //       goto TERMINATE_FREESTRINGS;
      //    }

      //    switch( option[1] )
      //    {
      //    case '-' :
      //       {
      //          option = &option[2];
      //          // This needs to be treated specially
      //          // --loadset=<setfile> : load parameters from settings file
      //          else if( strncmp(option, "loadset=", 8) == 0 )
      //          {
      //             if( loadsetname == nullptr )
      //             {
      //                char* filename = &option[8];
      //                loadsetname = new char[strlen(filename) + 1];
      //                spxSnprintf(loadsetname, strlen(filename) + 1, "%s", filename);
      //                if( !soplex->loadSettingsFile(loadsetname) )
      //                {
      //                   printUsage(argv, optidx);
      //                   returnValue = 1;
      //                   goto TERMINATE_FREESTRINGS;
      //                }
      //                else
      //                {
      //                   // we need to start parsing again because some command line parameters might have been overwritten
      //                   optidx = 0;
      //                }
      //             }
      //          }
      //          // --readmode=<value> : choose reading mode for <lpfile> (0* - floating-point, 1 - rational)
      //          else if( strncmp(option, "readmode=", 9) == 0 )
      //          {
      //             if( !soplex->setIntParam(soplex->READMODE, option[9] - '0') )
      //             {
      //                printUsage(argv, optidx);
      //                returnValue = 1;
      //                goto TERMINATE_FREESTRINGS;
      //             }
      //          }
      //          // --solvemode=<value> : choose solving mode (0* - floating-point solve, 1 - auto, 2 - force iterative refinement)
      //          else if( strncmp(option, "solvemode=", 10) == 0 )
      //          {
      //             if( !soplex->setIntParam(soplex->SOLVEMODE, option[10] - '0') )
      //             {
      //                printUsage(argv, optidx);
      //                returnValue = 1;
      //                goto TERMINATE_FREESTRINGS;
      //             }
      //             // if the LP is parsed rationally and might be solved rationally, we choose automatic syncmode such that
      //             // the rational LP is kept after reading
      //             else if( soplex->intParam(soplex->READMODE) == soplex->READMODE_RATIONAL
      //                && soplex->intParam(soplex->SOLVEMODE) != soplex->SOLVEMODE_REAL )
      //             {
      //                soplex->setIntParam(soplex->SYNCMODE, soplex->SYNCMODE_AUTO);
      //             }
      //          }
      //          // --extsol=<value> : external solution for soplex to use for validation
      //          else if( strncmp(option, "extsol=", 7) == 0 )
      //          {
      //             char* input = &option[7];
      //              if( !validation->updateExternalSolution(input) )
      //             {
      //                printUsage(argv, optidx);
      //                returnValue = 1;
      //                goto TERMINATE_FREESTRINGS;
      //             }
      //          }
      //          // --<type>:<name>=<val> :  change parameter value using syntax of settings file entries
      //          else if( !soplex->parseSettingsString(option) )
      //          {
      //             printUsage(argv, optidx);
      //             returnValue = 1;
      //             goto TERMINATE_FREESTRINGS;
      //          }
      //          break;
      //       }

      //    case 't' :
      //       // -t<s> : set time limit to <s> seconds
      //      if( !soplex->setRealParam(soplex->TIMELIMIT, atoi(&option[2])) )
      //       {
      //          printUsage(argv, optidx);
      //          returnValue = 1;
      //          goto TERMINATE_FREESTRINGS;
      //       }
      //       break;

      //    case 'i' :
      //       // -i<n> : set iteration limit to <n>
      //      if( !soplex->setIntParam(soplex->ITERLIMIT, atoi(&option[2])) )
      //       {
      //          printUsage(argv, optidx);
      //          returnValue = 1;
      //          goto TERMINATE_FREESTRINGS;
      //       }
      //       break;

      //    case 'f' :
      //       // -f<eps> : set primal feasibility tolerance to <eps>
      //      if( !soplex->setRealParam(soplex->FEASTOL, atof(&option[2])) )
      //       {
      //          printUsage(argv, optidx);
      //          returnValue = 1;
      //          goto TERMINATE_FREESTRINGS;
      //       }
      //       break;

      //    case 'o' :
      //       // -o<eps> : set dual feasibility (optimality) tolerance to <eps>
      //      if( !soplex->setRealParam(soplex->OPTTOL, atof(&option[2])) )
      //       {
      //          printUsage(argv, optidx);
      //          returnValue = 1;
      //          goto TERMINATE_FREESTRINGS;
      //       }
      //       break;

      //    case 'l' :
      //       // l<eps> : set validation tolerance to <eps>
      //       if( !validation->updateValidationTolerance(&option[2]) )
      //       {
      //          printUsage(argv, optidx);
      //          returnValue = 1;
      //          goto TERMINATE_FREESTRINGS;
      //       }
      //       break;

      //    case 's' :
      //       // -s<value> : choose simplifier/presolver (0 - off, 1* - auto)
      //      if( !soplex->setIntParam(soplex->SIMPLIFIER, option[2] - '0') )
      //       {
      //          printUsage(argv, optidx);
      //          returnValue = 1;
      //          goto TERMINATE_FREESTRINGS;
      //       }
      //       break;

      //    case 'g' :
      //       // -g<value> : choose scaling (0 - off, 1 - uni-equilibrium, 2* - bi-equilibrium, 3 - geometric, 4 - iterated geometric,  5 - least squares, 6 - geometric-equilibrium)
      //       if( !soplex->setIntParam(soplex->SCALER, option[2] - '0') )
      //       {
      //          printUsage(argv, optidx);
      //          returnValue = 1;
      //          goto TERMINATE_FREESTRINGS;
      //       }
      //       break;

      //    case 'p' :
      //       // -p<value> : choose pricing (0* - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - steep)
      //       if( !soplex->setIntParam(soplex->PRICER, option[2] - '0') )
      //       {
      //          printUsage(argv, optidx);
      //          returnValue = 1;
      //          goto TERMINATE_FREESTRINGS;
      //       }
      //       break;

      //    case 'r' :
      //       // -r<value> : choose ratio tester (0 - textbook, 1 - harris, 2* - fast, 3 - boundflipping)
      //       if( !soplex->setIntParam(soplex->RATIOTESTER, option[2] - '0') )
      //       {
      //          printUsage(argv, optidx);
      //          returnValue = 1;
      //          goto TERMINATE_FREESTRINGS;
      //       }
      //       break;

      //    case 'v' :
      //       // -v<level> : set verbosity to <level> (0 - error, 3 - normal, 5 - high)
      //       if( !soplex->setIntParam(soplex->VERBOSITY, option[2] - '0') )
      //       {
      //          printUsage(argv, optidx);
      //          returnValue = 1;
      //          goto TERMINATE_FREESTRINGS;
      //       }
      //       break;

      //    case 'x' :
      //       // -x : print primal solution
      //       printPrimal = true;
      //       break;

      //    case 'X' :
      //       // -X : print primal solution with rationals
      //       printPrimalRational = true;
      //       break;

      //    case 'y' :
      //       // -y : print dual multipliers
      //       printDual = true;
      //       break;

      //    case 'Y' :
      //       // -Y : print dual multipliers with rationals
      //       printDualRational = true;
      //       break;

      //    case 'q' :
      //       // -q : display detailed statistics
      //       displayStatistics = true;
      //       break;

      //    case 'c' :
      //       // -c : perform final check of optimal solution in original problem
      //       checkSol = true;
      //       break;

      //    case 'h' :
      //       // -h : display all parameters
      //       if( !soplex->saveSettingsFile(0, false) )
      //       {
      //          MSG_ERROR( std::cerr << "Error printing parameters\n" );
      //       }
      //       break;

      //       //lint -fallthrough
      //    default :
      //       {
      //          printUsage(argv, optidx);
      //          returnValue = 1;
      //          goto TERMINATE_FREESTRINGS;
      //       }
      //    }
      // }

      MSG_INFO1( soplex->spxout, soplex->printUserSettings(); )

        // TODO: How is the following code supposed to work?
      // no LP file was given and no settings files are written
        if( lpfilename.empty() && savesetname.empty() && diffsetname.empty())
      {
        // TODO: The printUsage should be different, also the GOTO macro
         printUsage(argv, 0);
         returnValue = 1;
         goto TERMINATE_FREESTRINGS;
      }

      // ensure that syncmode is not manual
      if( soplex->intParam(soplex->SYNCMODE) == soplex->SYNCMODE_MANUAL )
      {
         MSG_ERROR( std::cerr << "Error: manual synchronization is invalid on command line.  Change parameter int:syncmode.\n" );
         returnValue = 1;
         goto TERMINATE_FREESTRINGS;
      }

      // save settings files
      if(!savesetname.empty())
      {
         MSG_INFO1( soplex->spxout, soplex->spxout << "Saving parameters to settings file <" << savesetname << "> . . .\n" );
         if( !soplex->saveSettingsFile(savesetname, false) )
         {
            MSG_ERROR( std::cerr << "Error writing parameters to file <" << savesetname << ">\n" );
         }
      }
      if(!diffsetname.empty())
      {
         MSG_INFO1( soplex->spxout, soplex->spxout << "Saving modified parameters to settings file <" << diffsetname << "> . . .\n" );
         if( !soplex->saveSettingsFile(diffsetname, true) )
         {
            MSG_ERROR( std::cerr << "Error writing modified parameters to file <" << diffsetname << ">\n" );
         }
      }

      // no LP file given: exit after saving settings
      if(!lpfilename.empty())
      {
        if(!loadsetname.empty() || !savesetname.empty() || !diffsetname.empty())
         {
            MSG_INFO1( soplex->spxout, soplex->spxout << "\n" );
         }
         goto TERMINATE_FREESTRINGS;
      }

      // measure time for reading LP file and basis file
      readingTime->start();

      // if the LP is parsed rationally and might be solved rationally, we choose automatic syncmode such that
      // the rational LP is kept after reading
      if( soplex->intParam(soplex->READMODE) == soplex->READMODE_RATIONAL
         && soplex->intParam(soplex->SOLVEMODE) != soplex->SOLVEMODE_REAL )
      {
         soplex->setIntParam(soplex->SYNCMODE, soplex->SYNCMODE_AUTO);
      }

      // read LP from input file
      MSG_INFO1( soplex->spxout, soplex->spxout << "Reading "
         << (soplex->intParam(soplex->READMODE) == soplex->READMODE_REAL ? "(real)" : "(rational)")
         << " LP file <" << lpfilename << "> . . .\n" );

      if( !soplex->readFile(lpfilename, &rownames, &colnames) )
      {
         MSG_ERROR( std::cerr << "Error while reading file <" << lpfilename << ">.\n" );
         returnValue = 1;
         goto TERMINATE_FREESTRINGS;
      }

      // write LP if specified
      if(!writefilename.empty())
      {
        if( !soplex->writeFile(writefilename.c_str(), &rownames, &colnames) )
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

      // write dual LP if specified
      if(!writedualfilename.empty())
      {
        if( !soplex->writeDualFileReal(writedualfilename.c_str(), &rownames, &colnames) )
         {
            MSG_ERROR( std::cerr << "Error while writing dual file <" << writedualfilename << ">.\n\n" );
            returnValue = 1;
            goto TERMINATE_FREESTRINGS;
         }
         else
         {
            MSG_INFO1( soplex->spxout, soplex->spxout << "Written dual LP to file <" << writedualfilename << ">.\n\n" );
         }
      }

      // read basis file if specified
      if(!readbasname.empty())
      {
         MSG_INFO1( soplex->spxout, soplex->spxout << "Reading basis file <" << readbasname << "> . . . " );
         if( !soplex->readBasisFile(readbasname.c_str(), &rownames, &colnames) )
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

      MSG_INFO1( soplex->spxout, soplex->spxout << "LP has " << soplex->numRows() << " rows "
                 << soplex->numCols() << " columns and " << soplex->numNonzeros() << " nonzeros.\n\n" );

      // solve the LP
      soplex->optimize();

      // print solution, check solution, and display statistics
      printPrimalSolution(*soplex, colnames, rownames, printPrimal, printPrimalRational);
      printDualSolution(*soplex, colnames, rownames, printDual, printDualRational);

      if( checkSol )
        checkSolution<R>(*soplex); // The type needs to get fixed here

      if( displayStatistics )
      {
         MSG_INFO1( soplex->spxout, soplex->spxout << "Statistics\n==========\n\n" );
         soplex->printStatistics(soplex->spxout.getStream(SPxOut::INFO1));
      }

      if(validation->validate)
         validation->validateSolveReal(*soplex);

      // write basis file if specified
      if(!writebasname.empty())
      {
         if( !soplex->hasBasis() )
         {
            MSG_WARNING( soplex->spxout, soplex->spxout << "No basis information available.  Could not write file <" << writebasname << ">\n\n" );
         }
         else if( !soplex->writeBasisFile(writebasname.c_str(), &rownames, &colnames) )
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
   catch( const SPxException& x ) // There could be an exception from boost?
   {
      MSG_ERROR( std::cerr << "Exception caught: " << x.what() << "\n" );
      returnValue = 1;
      goto TERMINATE_FREESTRINGS;
   }

 TERMINATE:
   // because EGlpNumClear() calls mpq_clear() for all mpq_t variables, we need to destroy all objects of class Rational
   // beforehand; hence all Rational objects and all data that uses Rational objects must be allocated dynamically via
   // spx_alloc() and freed here; disabling the list memory is crucial
 if( nullptr != soplex )
 {
    soplex->~SoPlexBase();
    spx_free(soplex);
 }
 if( nullptr != validation )
 {
    validation->~Validation();
    spx_free(validation);
 }
 Rational::disableListMem();
 EGlpNumClear();
 if( nullptr != readingTime )
 {
    readingTime->~Timer();
    spx_free(readingTime);
 }

   return returnValue;
}



}
