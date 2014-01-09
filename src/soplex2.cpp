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

#include <iostream>
#include <assert.h>

#include "soplex2.h"
#include "spxfileio.h"
#include "statistics.h"
#include "mpsinput.h"
#include "limits.h"

/// maximum length of lines in settings file
#define SET_MAX_LINE_LEN 500

#ifdef _MSC_VER
#define strncasecmp strnicmp
#endif

namespace soplex
{
   /// class of parameter settings
   class SoPlex2::Settings
   {
   public:
      /// array of names for boolean parameters
      static std::string _boolParamName[SoPlex2::BOOLPARAM_COUNT];

      /// array of names for integer parameters
      static std::string _intParamName[SoPlex2::INTPARAM_COUNT];

      /// array of names for real parameters
      static std::string _realParamName[SoPlex2::REALPARAM_COUNT];

      /// array of names for rational parameters
      static std::string _rationalParamName[SoPlex2::RATIONALPARAM_COUNT];

      /// array of descriptions for boolean parameters
      static std::string _boolParamDescription[SoPlex2::BOOLPARAM_COUNT];

      /// array of descriptions for integer parameters
      static std::string _intParamDescription[SoPlex2::INTPARAM_COUNT];

      /// array of descriptions for real parameters
      static std::string _realParamDescription[SoPlex2::REALPARAM_COUNT];

      /// array of descriptions for rational parameters
      static std::string _rationalParamDescription[SoPlex2::RATIONALPARAM_COUNT];

      /// array of default values for boolean parameters
      static bool _boolParamDefault[SoPlex2::BOOLPARAM_COUNT];

      /// array of default values for integer parameters
      static int _intParamDefault[SoPlex2::INTPARAM_COUNT];

      /// array of default values for real parameters
      static Real _realParamDefault[SoPlex2::REALPARAM_COUNT];

      /// array of default values for rational parameters
      static Rational _rationalParamDefault[SoPlex2::RATIONALPARAM_COUNT];

      /// array of lower bounds for real parameter values
      static Real _realParamLower[SoPlex2::REALPARAM_COUNT];

      /// array of upper bounds for real parameter values
      static Real _realParamUpper[SoPlex2::REALPARAM_COUNT];

      /// array of lower bounds for rational parameter values
      static Rational _rationalParamLower[SoPlex2::RATIONALPARAM_COUNT];

      /// array of upper bounds for rational parameter values
      static Rational _rationalParamUpper[SoPlex2::RATIONALPARAM_COUNT];

      /// have static arrays been initialized?
      static bool _defaultsAndBoundsInitialized;

      /// array of current boolean parameter values
      bool _boolParamValues[SoPlex2::BOOLPARAM_COUNT];

      /// array of current integer parameter values
      int _intParamValues[SoPlex2::INTPARAM_COUNT];

      /// array of current real parameter values
      Real _realParamValues[SoPlex2::REALPARAM_COUNT];

      /// array of current rational parameter values
      Rational _rationalParamValues[SoPlex2::RATIONALPARAM_COUNT];

      /// default constructor initializing default settings
      Settings()
      {
         if( !_defaultsAndBoundsInitialized )
         {
            // should partial pricing be used?
            _boolParamName[SoPlex2::PARTIAL_PRICING] = "partial_pricing";
            _boolParamDescription[SoPlex2::PARTIAL_PRICING] = "should partial pricing be used?";
            _boolParamDefault[SoPlex2::PARTIAL_PRICING] = false;

            // should lifting be used to reduce range of nonzero matrix coefficients?
            _boolParamName[SoPlex2::LIFTING] = "lifting";
            _boolParamDescription[SoPlex2::LIFTING] = "should lifting be used to reduce range of nonzero matrix coefficients?";
            _boolParamDefault[SoPlex2::LIFTING] = false;

            // objective sense
            _intParamName[SoPlex2::OBJSENSE] = "objsense";
            _intParamDescription[SoPlex2::OBJSENSE] = "objective sense (-1 - minimize, +1 - maximize)";
            _intParamDefault[SoPlex2::OBJSENSE] = SoPlex2::OBJSENSE_MAXIMIZE;

            // type of computational form, i.e., column or row representation
            _intParamName[SoPlex2::REPRESENTATION] = "representation";
            _intParamDescription[SoPlex2::REPRESENTATION] = "type of computational form (0 - auto, 1 - column representation, 2 - row representation)";
            _intParamDefault[SoPlex2::REPRESENTATION] = SoPlex2::REPRESENTATION_AUTO;

            // type of algorithm, i.e., enter or leave
            _intParamName[SoPlex2::ALGORITHM] = "algorithm";
            _intParamDescription[SoPlex2::ALGORITHM] = "type of algorithm (0 - enter, 1 - leave)";
            _intParamDefault[SoPlex2::ALGORITHM] = SoPlex2::ALGORITHM_LEAVE;

            // type of LU update
            _intParamName[SoPlex2::FACTOR_UPDATE_TYPE] = "factor_update_type";
            _intParamDescription[SoPlex2::FACTOR_UPDATE_TYPE] = "type of LU update (0 - eta update, 1 - Forrest-Tomlin update)";
            _intParamDefault[SoPlex2::FACTOR_UPDATE_TYPE] = SoPlex2::FACTOR_UPDATE_TYPE_FT;

            ///@todo which value?
            // maximum number of updates without fresh factorization
            _intParamName[SoPlex2::FACTOR_UPDATE_MAX] = "factor_update_max";
            _intParamDescription[SoPlex2::FACTOR_UPDATE_MAX] = "maximum number of LU updates without fresh factorization";
            _intParamDefault[SoPlex2::FACTOR_UPDATE_MAX] = 200;

            // iteration limit (-1 if unlimited)
            _intParamName[SoPlex2::ITERLIMIT] = "iterlimit";
            _intParamDescription[SoPlex2::ITERLIMIT] = "iteration limit (-1 - no limit)";
            _intParamDefault[SoPlex2::ITERLIMIT] = -1;

            // refinement limit (-1 if unlimited)
            _intParamName[SoPlex2::REFLIMIT] = "reflimit";
            _intParamDescription[SoPlex2::REFLIMIT] = "refinement limit (-1 - no limit)";
            _intParamDefault[SoPlex2::REFLIMIT] = -1;

            // stalling refinement limit (-1 if unlimited)
            _intParamName[SoPlex2::STALLREFLIMIT] = "stallreflimit";
            _intParamDescription[SoPlex2::STALLREFLIMIT] = "stalling refinement limit (-1 - no limit)";
            _intParamDefault[SoPlex2::STALLREFLIMIT] = -1;

            // display frequency
            _intParamName[SoPlex2::DISPLAYFREQ] = "displayfreq";
            _intParamDescription[SoPlex2::DISPLAYFREQ] = "display frequency";
            _intParamDefault[SoPlex2::DISPLAYFREQ] = 100;

            // verbosity level
            _intParamName[SoPlex2::VERBOSITY] = "verbosity";
            _intParamDescription[SoPlex2::VERBOSITY] = "verbosity level (0 - error, 1 - warning, 2 - debug, 3 - normal, 4 - high, 5 - full)";
            _intParamDefault[SoPlex2::VERBOSITY] = SoPlex2::VERBOSITY_NORMAL;

            // type of simplifier
            _intParamName[SoPlex2::SIMPLIFIER] = "simplifier";
            _intParamDescription[SoPlex2::SIMPLIFIER] = "simplifier (0 - off, 1 - auto)";
            _intParamDefault[SoPlex2::SIMPLIFIER] = SoPlex2::SIMPLIFIER_AUTO;

            // type of scaler
            _intParamName[SoPlex2::SCALER] = "scaler";
            _intParamDescription[SoPlex2::SCALER] = "scaling (0 - off, 1 - uni-equilibrium, 2 - bi-equilibrium, 3 - geometric, 4 - iterated geometric)";
            _intParamDefault[SoPlex2::SCALER] = SoPlex2::SCALER_BIEQUI;

            // type of starter used to create crash basis
            _intParamName[SoPlex2::STARTER] = "starter";
            _intParamDescription[SoPlex2::STARTER] = "crash basis generated when starting from scratch (0 - none, 1 - weight, 2 - sum, 3 - vector)";
            _intParamDefault[SoPlex2::STARTER] = SoPlex2::STARTER_OFF;

            // type of pricer
            _intParamName[SoPlex2::PRICER] = "pricer";
            _intParamDescription[SoPlex2::PRICER] = "pricing method (0 - auto, 1 - dantzig, 2 - parmult, 3 - devex, 4 - quicksteep, 5 - hybrid quicksteep/parmult)";
            _intParamDefault[SoPlex2::PRICER] = SoPlex2::PRICER_QUICKSTEEP;

            // type of ratio test
            _intParamName[SoPlex2::RATIOTESTER] = "ratiotester";
            _intParamDescription[SoPlex2::RATIOTESTER] = "method for ratio test (0 - textbook, 1 - harris, 2 - fast, 3 - boundflipping)";
            _intParamDefault[SoPlex2::RATIOTESTER] = SoPlex2::RATIOTESTER_FAST;

            // mode for synchronizing real and rational LP
            _intParamName[SoPlex2::SYNCMODE] = "syncmode";
            _intParamDescription[SoPlex2::SYNCMODE] = "mode for synchronizing real and rational LP (0 - store only real LP, 1 - auto, 2 - manual)";
            _intParamDefault[SoPlex2::SYNCMODE] = SoPlex2::SYNCMODE_ONLYREAL;

            // mode for reading LP files
            _intParamName[SoPlex2::READMODE] = "readmode";
            _intParamDescription[SoPlex2::READMODE] = "mode for reading LP files (0 - floating-point, 1 - rational)";
            _intParamDefault[SoPlex2::READMODE] = SoPlex2::READMODE_REAL;

            // mode for iterative refinement strategy
            _intParamName[SoPlex2::SOLVEMODE] = "solvemode";
            _intParamDescription[SoPlex2::SOLVEMODE] = "mode for iterative refinement strategy (0 - floating-point solve, 1 - auto, 2 - force iterative refinement)";
            _intParamDefault[SoPlex2::SOLVEMODE] = SoPlex2::SOLVEMODE_REAL;

            ///@todo define suitable values depending on Real type
            // general zero tolerance
            _realParamName[SoPlex2::EPSILON_ZERO] = "epsilon_zero";
            _realParamDescription[SoPlex2::EPSILON_ZERO] = "general zero tolerance";
            _realParamLower[SoPlex2::EPSILON_ZERO] = DEFAULT_EPS_ZERO;
            _realParamUpper[SoPlex2::EPSILON_ZERO] = DEFAULT_EPS_ZERO;
            _realParamDefault[SoPlex2::EPSILON_ZERO] = DEFAULT_EPS_ZERO;

            ///@todo define suitable values depending on Real type
            // zero tolerance used in factorization
            _realParamName[SoPlex2::EPSILON_FACTORIZATION] = "epsilon_factorization";
            _realParamDescription[SoPlex2::EPSILON_FACTORIZATION] = "zero tolerance used in factorization";
            _realParamLower[SoPlex2::EPSILON_FACTORIZATION] = DEFAULT_EPS_FACTOR;
            _realParamUpper[SoPlex2::EPSILON_FACTORIZATION] = DEFAULT_EPS_FACTOR;
            _realParamDefault[SoPlex2::EPSILON_FACTORIZATION] = DEFAULT_EPS_FACTOR;

            ///@todo define suitable values depending on Real type
            // zero tolerance used in update of the factorization
            _realParamName[SoPlex2::EPSILON_UPDATE] = "epsilon_update";
            _realParamDescription[SoPlex2::EPSILON_UPDATE] = "zero tolerance used in update of the factorization";
            _realParamLower[SoPlex2::EPSILON_UPDATE] = DEFAULT_EPS_UPDATE;
            _realParamUpper[SoPlex2::EPSILON_UPDATE] = DEFAULT_EPS_UPDATE;
            _realParamDefault[SoPlex2::EPSILON_UPDATE] = DEFAULT_EPS_UPDATE;

            ///@todo define suitable values depending on Real type
            // infinity threshold
            _realParamName[SoPlex2::INFTY] = "infty";
            _realParamDescription[SoPlex2::INFTY] = "infinity threshold";
            _realParamLower[SoPlex2::INFTY] = DEFAULT_INFINITY;
            _realParamUpper[SoPlex2::INFTY] = DEFAULT_INFINITY;
            _realParamDefault[SoPlex2::INFTY] = DEFAULT_INFINITY;

            // time limit in seconds (INFTY if unlimited)
            _realParamName[SoPlex2::TIMELIMIT] = "timelimit";
            _realParamDescription[SoPlex2::TIMELIMIT] = "time limit in seconds";
            _realParamLower[SoPlex2::TIMELIMIT] = 0.0;
            _realParamUpper[SoPlex2::TIMELIMIT] = DEFAULT_INFINITY;
            _realParamDefault[SoPlex2::TIMELIMIT] = DEFAULT_INFINITY;

            // lower limit on objective value
            _realParamName[SoPlex2::OBJLIMIT_LOWER] = "objlimit_lower";
            _realParamDescription[SoPlex2::OBJLIMIT_LOWER] = "lower limit on objective value";
            _realParamLower[SoPlex2::OBJLIMIT_LOWER] = -_realParamLower[SoPlex2::INFTY];
            _realParamUpper[SoPlex2::OBJLIMIT_LOWER] = _realParamLower[SoPlex2::INFTY];
            _realParamDefault[SoPlex2::OBJLIMIT_LOWER] = -_realParamLower[SoPlex2::INFTY];

            // upper limit on objective value
            _realParamName[SoPlex2::OBJLIMIT_UPPER] = "objlimit_upper";
            _realParamDescription[SoPlex2::OBJLIMIT_UPPER] = "upper limit on objective value";
            _realParamLower[SoPlex2::OBJLIMIT_UPPER] = -_realParamLower[SoPlex2::INFTY];
            _realParamUpper[SoPlex2::OBJLIMIT_UPPER] = _realParamLower[SoPlex2::INFTY];
            _realParamDefault[SoPlex2::OBJLIMIT_UPPER] = _realParamLower[SoPlex2::INFTY];

            // working tolerance for feasibility in floating-point solver during iterative refinement
            _realParamName[SoPlex2::FPFEASTOL] = "fpfeastol";
            _realParamDescription[SoPlex2::FPFEASTOL] = "working tolerance for feasibility in floating-point solver during iterative refinement";
            _realParamLower[SoPlex2::FPFEASTOL] = 1e-12;
            _realParamUpper[SoPlex2::FPFEASTOL] = 1.0;
            _realParamDefault[SoPlex2::FPFEASTOL] = 1e-6;

            // working tolerance for optimality in floating-point solver during iterative refinement
            _realParamName[SoPlex2::FPOPTTOL] = "fpopttol";
            _realParamDescription[SoPlex2::FPOPTTOL] = "working tolerance for optimality in floating-point solver during iterative refinement";
            _realParamLower[SoPlex2::FPOPTTOL] = 1e-12;
            _realParamUpper[SoPlex2::FPOPTTOL] = 1.0;
            _realParamDefault[SoPlex2::FPOPTTOL] = 1e-6;

            // maximum increase of scaling factors between refinements
            _realParamName[SoPlex2::MAXSCALEINCR] = "maxscaleincr";
            _realParamDescription[SoPlex2::MAXSCALEINCR] = "maximum increase of scaling factors between refinements";
            _realParamLower[SoPlex2::MAXSCALEINCR] = 1.0;
            _realParamUpper[SoPlex2::MAXSCALEINCR] = DEFAULT_INFINITY;
            _realParamDefault[SoPlex2::MAXSCALEINCR] = 1e12;

            // lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)
            _realParamName[SoPlex2::LIFTMINVAL] = "liftminval";
            _realParamDescription[SoPlex2::LIFTMINVAL] = "lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)";
            _realParamLower[SoPlex2::LIFTMINVAL] = 0.0;
            _realParamUpper[SoPlex2::LIFTMINVAL] = 0.1;
            _realParamDefault[SoPlex2::LIFTMINVAL] = 0.000976562; // = 1/1024

            // upper threshold in lifting (nonzero matrix coefficients with larger absolute value will be reformulated)
            _realParamName[SoPlex2::LIFTMAXVAL] = "liftmaxval";
            _realParamDescription[SoPlex2::LIFTMAXVAL] = "lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)";
            _realParamLower[SoPlex2::LIFTMAXVAL] = 10.0;
            _realParamUpper[SoPlex2::LIFTMAXVAL] = DEFAULT_INFINITY;
            _realParamDefault[SoPlex2::LIFTMAXVAL] = 1024.0;

            // primal feasibility tolerance
            _rationalParamName[SoPlex2::FEASTOL] = "feastol";
            _rationalParamDescription[SoPlex2::FEASTOL] = "primal feasibility tolerance";
            _rationalParamLower[SoPlex2::FEASTOL] = 0.0;
            _rationalParamUpper[SoPlex2::FEASTOL] = 1.0;
            _rationalParamDefault[SoPlex2::FEASTOL] = Rational(1)/Rational(1000000);

            // dual feasibility tolerance
            _rationalParamName[SoPlex2::OPTTOL] = "opttol";
            _rationalParamDescription[SoPlex2::OPTTOL] = "dual feasibility tolerance";
            _rationalParamLower[SoPlex2::OPTTOL] = 0.0;
            _rationalParamUpper[SoPlex2::OPTTOL] = 1.0;
            _rationalParamDefault[SoPlex2::OPTTOL] = Rational(1)/Rational(1000000);

            _defaultsAndBoundsInitialized = true;
         }

         for( int i = 0; i < SoPlex2::BOOLPARAM_COUNT; i++ )
            _boolParamValues[i] = _boolParamDefault[i];

         for( int i = 0; i < SoPlex2::INTPARAM_COUNT; i++ )
            _intParamValues[i] = _intParamDefault[i];

         for( int i = 0; i < SoPlex2::REALPARAM_COUNT; i++ )
            _realParamValues[i] = _realParamDefault[i];

         for( int i = 0; i < SoPlex2::RATIONALPARAM_COUNT; i++ )
            _rationalParamValues[i] = _rationalParamDefault[i];
      }

      /// copy constructor
      Settings(const Settings& settings)
      {
         *this = settings;
      }

      /// assignment operator
      Settings& operator=(const Settings& settings)
      {
         for( int i = 0; i < SoPlex2::BOOLPARAM_COUNT; i++ )
            _boolParamValues[i] = settings._boolParamValues[i];

         for( int i = 0; i < SoPlex2::INTPARAM_COUNT; i++ )
            _intParamValues[i] = settings._intParamValues[i];

         for( int i = 0; i < SoPlex2::REALPARAM_COUNT; i++ )
            _realParamValues[i] = settings._realParamValues[i];

         for( int i = 0; i < SoPlex2::RATIONALPARAM_COUNT; i++ )
            _rationalParamValues[i] = settings._rationalParamValues[i];

         return *this;
      }
   };



   bool SoPlex2::Settings::_defaultsAndBoundsInitialized = false;



   std::string SoPlex2::Settings::_boolParamName[SoPlex2::BOOLPARAM_COUNT];
   std::string SoPlex2::Settings::_boolParamDescription[SoPlex2::BOOLPARAM_COUNT];
   bool SoPlex2::Settings::_boolParamDefault[SoPlex2::BOOLPARAM_COUNT];



   std::string SoPlex2::Settings::_intParamName[SoPlex2::INTPARAM_COUNT];
   std::string SoPlex2::Settings::_intParamDescription[SoPlex2::INTPARAM_COUNT];
   int SoPlex2::Settings::_intParamDefault[SoPlex2::INTPARAM_COUNT];



   std::string SoPlex2::Settings::_realParamName[SoPlex2::REALPARAM_COUNT];
   std::string SoPlex2::Settings::_realParamDescription[SoPlex2::REALPARAM_COUNT];
   Real SoPlex2::Settings::_realParamLower[SoPlex2::REALPARAM_COUNT];
   Real SoPlex2::Settings::_realParamUpper[SoPlex2::REALPARAM_COUNT];
   Real SoPlex2::Settings::_realParamDefault[SoPlex2::REALPARAM_COUNT];



   std::string SoPlex2::Settings::_rationalParamName[SoPlex2::RATIONALPARAM_COUNT];
   std::string SoPlex2::Settings::_rationalParamDescription[SoPlex2::RATIONALPARAM_COUNT];
   Rational SoPlex2::Settings::_rationalParamLower[SoPlex2::RATIONALPARAM_COUNT];
   Rational SoPlex2::Settings::_rationalParamUpper[SoPlex2::RATIONALPARAM_COUNT];
   Rational SoPlex2::Settings::_rationalParamDefault[SoPlex2::RATIONALPARAM_COUNT];



   /// default constructor
   SoPlex2::SoPlex2()
      : _rationalLP(0)
      , _statistics(0)
      , _currentSettings(0)
      , _scalerUniequi(false)
      , _scalerBiequi(true)
      , _scalerGeo1(1)
      , _scalerGeo8(8)
      , _simplifier(0)
      , _scaler(0)
      , _starter(0)
      , _status(SPxSolver::UNKNOWN)
      , _hasBasis(false)
      , _hasSolReal(false)
      , _hasSolRational(false)
   {
      // give lu factorization to solver
      _solver.setSolver(&_slufactor);

      // the real LP is initially stored in the solver; the rational LP is constructed, when the parameter SYNCMODE is
      // initialized in setSettings() below
      _realLP = &_solver;
      _isRealLPLoaded = true;

      // initialize statistics
      spx_alloc(_statistics);
      _statistics = new (_statistics) Statistics();

      // initialize parameter settings to default
      spx_alloc(_currentSettings);
      _currentSettings = new (_currentSettings) Settings();
      setSettings(*_currentSettings, true, true);

      assert(_isConsistent());
   }



   /// assignment operator
   SoPlex2& SoPlex2::operator=(const SoPlex2& rhs)
   {
      if( this != &rhs )
      {
         // copy statistics
         *_statistics = *(rhs._statistics);

         // copy settings
         *_currentSettings = *(rhs._currentSettings);

         // copy solver components
         _solver = rhs._solver;
         _slufactor = rhs._slufactor;
         _simplifierMainSM = rhs._simplifierMainSM;
         _scalerUniequi = rhs._scalerUniequi;
         _scalerBiequi = rhs._scalerBiequi;
         _scalerGeo1 = rhs._scalerGeo1;
         _scalerGeo8 = rhs._scalerGeo8;
         _starterWeight = rhs._starterWeight;
         _starterSum = rhs._starterSum;
         _starterVector = rhs._starterVector;
         _pricerDantzig = rhs._pricerDantzig;
         _pricerParMult = rhs._pricerParMult;
         _pricerDevex = rhs._pricerDevex;
         _pricerQuickSteep = rhs._pricerQuickSteep;
         _pricerSteep = rhs._pricerSteep;
         _pricerHybrid = rhs._pricerHybrid;
         _ratiotesterTextbook = rhs._ratiotesterTextbook;
         _ratiotesterHarris = rhs._ratiotesterHarris;
         _ratiotesterFast = rhs._ratiotesterFast;
         _ratiotesterBoundFlipping = rhs._ratiotesterBoundFlipping;

         // copy solution data
         _status = rhs._status;
         _basisStatusRows = rhs._basisStatusRows;
         _basisStatusCols = rhs._basisStatusCols;

         if( rhs._hasSolReal )
            _solReal = rhs._solReal;

         if( rhs._hasSolRational )
            _solRational = rhs._solRational;

         // initialize pointers for simplifier, scaler, and starter
         setIntParam(SoPlex2::SIMPLIFIER, intParam(SoPlex2::SIMPLIFIER), true, true);
         setIntParam(SoPlex2::SCALER, intParam(SoPlex2::SCALER), true, true);
         setIntParam(SoPlex2::STARTER, intParam(SoPlex2::STARTER), true, true);

         // copy real LP if different from the LP in the solver
         if( rhs._realLP != &(rhs._solver) )
         {
            _realLP = 0;
            spx_alloc(_realLP);
            _realLP = new (_realLP) SPxLPReal(*(rhs._realLP));
         }
         else
            _realLP = &_solver;

         // copy rational LP
         if( rhs._rationalLP == 0 )
         {
            assert(intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL);
            _rationalLP = 0;
         }
         else
         {
            assert(intParam(SoPlex2::SYNCMODE) != SYNCMODE_ONLYREAL);
            _rationalLP = 0;
            spx_alloc(_rationalLP);
            _rationalLP = new (_rationalLP) SPxLPRational(*rhs._rationalLP);
         }

         // copy boolean flags
         _isRealLPLoaded = rhs._isRealLPLoaded;
         _hasSolReal = rhs._hasSolReal;
         _hasSolRational = rhs._hasSolRational;
         _hasBasis = rhs._hasBasis;
      }

      assert(_isConsistent());

      return *this;
   }



   /// copy constructor
   SoPlex2::SoPlex2(const SoPlex2& rhs)
   {
      ///@todo improve performance by implementing a separate copy constructor
      *this = rhs;
   }



   /// destructor
   SoPlex2::~SoPlex2()
   {
      assert(_isConsistent());

      // free settings
      _currentSettings->~Settings();
      spx_free(_currentSettings);

      // free statistics
      _statistics->~Statistics();
      spx_free(_statistics);

      // free real LP if different from the LP in the solver
      assert(_realLP != 0);
      if( _realLP != &_solver )
      {
         _realLP->~SPxLPReal();
         spx_free(_realLP);
      }

      // free rational LP
      if( _rationalLP != 0 )
      {
         _rationalLP->~SPxLPRational();
         spx_free(_rationalLP);
      }
   }



   /// returns number of rows
   int SoPlex2::numRowsReal() const
   {
      assert(_realLP != 0);
      return _realLP->nRows();
   }



   /// returns number of columns
   int SoPlex2::numColsReal() const
   {
      assert(_realLP != 0);
      return _realLP->nCols();
   }



   /// returns number of nonzeros
   int SoPlex2::numNonzerosReal() const
   {
      assert(_realLP != 0);
      return _realLP->nNzos();
   }



   /// returns smallest non-zero element in absolute value
   Real SoPlex2::minAbsNonzeroReal() const
   {
      assert(_realLP != 0);
      return _realLP->minAbsNzo();
   }



   /// returns biggest non-zero element in absolute value
   Real SoPlex2::maxAbsNonzeroReal() const
   {
      assert(_realLP != 0);
      return _realLP->maxAbsNzo();
   }



   /// gets row \p i
   void SoPlex2::getRowReal(int i, LPRowReal& lprow) const
   {
      assert(_realLP != 0);
      _realLP->getRow(i, lprow);
   }



   /// gets rows \p start, ..., \p end.
   void SoPlex2::getRowsReal(int start, int end, LPRowSetReal& lprowset) const
   {
      assert(_realLP != 0);
      _realLP->getRows(start, end, lprowset);
   }



   /// returns vector of row \p i
   const SVectorReal& SoPlex2::rowVectorReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->rowVector(i);
   }



   /// returns right-hand side vector
   const VectorReal& SoPlex2::rhsReal() const
   {
      assert(_realLP != 0);
      return _realLP->rhs();
   }



   /// returns right-hand side of row \p i
   Real SoPlex2::rhsReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->rhs(i);
   }



   /// returns left-hand side vector
   const VectorReal& SoPlex2::lhsReal() const
   {
      assert(_realLP != 0);
      return _realLP->lhs();
   }



   /// returns left-hand side of row \p i
   Real SoPlex2::lhsReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->lhs(i);
   }



   /// returns inequality type of row \p i
   LPRowReal::Type SoPlex2::rowTypeReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->rowType(i);
   }



   /// gets column \p i
   void SoPlex2::getColReal(int i, LPColReal& lpcol) const
   {
      assert(_realLP != 0);
      return _realLP->getCol(i, lpcol);
   }



   /// gets columns \p start, ..., \p end
   void SoPlex2::getColsReal(int start, int end, LPColSetReal& lpcolset) const
   {
      assert(_realLP != 0);
      return _realLP->getCols(start, end, lpcolset);
   }



   /// returns vector of column \p i
   const SVectorReal& SoPlex2::colVectorReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->colVector(i);
   }



   /// returns upper bound vector
   const VectorReal& SoPlex2::upperReal() const
   {
      assert(_realLP != 0);
      return _realLP->upper();
   }



   /// returns upper bound of column \p i
   Real SoPlex2::upperReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->upper(i);
   }



   /// returns lower bound vector
   const VectorReal& SoPlex2::lowerReal() const
   {
      assert(_realLP != 0);
      return _realLP->lower();
   }



   /// returns lower bound of column \p i
   Real SoPlex2::lowerReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->lower(i);
   }



   /// gets objective function vector
   void SoPlex2::getObjReal(VectorReal& obj) const
   {
      assert(_realLP != 0);
      _realLP->getObj(obj);
   }



   /// returns objective value of column \p i
   Real SoPlex2::objReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->obj(i);
   }



   /// returns objective function vector after transformation to a maximization problem; since this is how it is stored
   /// internally, this is generally faster
   const VectorReal& SoPlex2::maxObjReal() const
   {
      assert(_realLP != 0);
      return _realLP->maxObj();
   }



   /// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
   /// stored internally, this is generally faster
   Real SoPlex2::maxObjReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->maxObj(i);
   }



   /// returns number of rows
   int SoPlex2::numRowsRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->nRows();
   }



   /// returns number of columns
   int SoPlex2::numColsRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->nCols();
   }



   /// returns number of nonzeros
   int SoPlex2::numNonzerosRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->nNzos();
   }



   /// returns smallest non-zero element in absolute value
   Rational SoPlex2::minAbsNonzeroRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->minAbsNzo();
   }



   /// returns biggest non-zero element in absolute value
   Rational SoPlex2::maxAbsNonzeroRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->maxAbsNzo();
   }



   /// gets row \p i
   void SoPlex2::getRowRational(int i, LPRowRational& lprow) const
   {
      assert(_rationalLP != 0);
      _rationalLP->getRow(i, lprow);
   }



   /// gets rows \p start, ..., \p end.
   void SoPlex2::getRowsRational(int start, int end, LPRowSetRational& lprowset) const
   {
      assert(_rationalLP != 0);
      _rationalLP->getRows(start, end, lprowset);
   }



   /// returns vector of row \p i
   const SVectorRational& SoPlex2::rowVectorRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->rowVector(i);
   }



   /// returns right-hand side vector
   const VectorRational& SoPlex2::rhsRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->rhs();
   }



   /// returns right-hand side of row \p i
   Rational SoPlex2::rhsRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->rhs(i);
   }



   /// returns left-hand side vector
   const VectorRational& SoPlex2::lhsRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->lhs();
   }



   /// returns left-hand side of row \p i
   Rational SoPlex2::lhsRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->lhs(i);
   }



   /// returns inequality type of row \p i
   LPRowRational::Type SoPlex2::rowTypeRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->rowType(i);
   }



   /// gets column \p i
   void SoPlex2::getColRational(int i, LPColRational& lpcol) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->getCol(i, lpcol);
   }



   /// gets columns \p start, ..., \p end
   void SoPlex2::getColsRational(int start, int end, LPColSetRational& lpcolset) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->getCols(start, end, lpcolset);
   }



   /// returns vector of column \p i
   const SVectorRational& SoPlex2::colVectorRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->colVector(i);
   }



   /// returns upper bound vector
   const VectorRational& SoPlex2::upperRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->upper();
   }



   /// returns upper bound of column \p i
   Rational SoPlex2::upperRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->upper(i);
   }



   /// returns lower bound vector
   const VectorRational& SoPlex2::lowerRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->lower();
   }



   /// returns lower bound of column \p i
   Rational SoPlex2::lowerRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->lower(i);
   }



   /// gets objective function vector
   void SoPlex2::getObjRational(VectorRational& obj) const
   {
      assert(_rationalLP != 0);
      _rationalLP->getObj(obj);
   }



   /// returns objective value of column \p i
   Rational SoPlex2::objRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->obj(i);
   }



   /// returns objective function vector after transformation to a maximization problem; since this is how it is stored
   /// internally, this is generally faster
   const VectorRational& SoPlex2::maxObjRational() const
   {
      assert(_rationalLP != 0);
      return _rationalLP->maxObj();
   }



   /// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
   /// stored internally, this is generally faster
   Rational SoPlex2::maxObjRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->maxObj(i);
   }



   /// adds a single row
   void SoPlex2::addRowReal(const LPRowReal& lprow)
   {
      assert(_realLP != 0);

      _addRowReal(lprow);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
        _rationalLP->addRow(lprow);

      _invalidateSolution();
   }



   /// adds multiple rows
   void SoPlex2::addRowsReal(const LPRowSetReal& lprowset)
   {
      assert(_realLP != 0);

      _addRowsReal(lprowset);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->addRows(lprowset);

      _invalidateSolution();
   }



   /// adds a single column
   void SoPlex2::addColReal(const LPColReal& lpcol)
   {
      assert(_realLP != 0);

      _addColReal(lpcol);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->addCol(lpcol);

      _invalidateSolution();
   }



   /// adds multiple columns
   void SoPlex2::addColsReal(const LPColSetReal& lpcolset)
   {
      assert(_realLP != 0);

      _addColsReal(lpcolset);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->addCols(lpcolset);

      _invalidateSolution();
   }



   /// replaces row \p i with \p lprow
   void SoPlex2::changeRowReal(int i, const LPRowReal& lprow)
   {
      assert(_realLP != 0);

      _changeRowReal(i, lprow);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeRow(i, lprow);

      _invalidateSolution();
   }



   /// changes left-hand side vector for constraints to \p lhs
   void SoPlex2::changeLhsReal(const VectorReal& lhs)
   {
      assert(_realLP != 0);

      _changeLhsReal(lhs);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeLhs(DVectorRational(lhs));

      _invalidateSolution();
   }



   /// changes left-hand side of row \p i to \p lhs
   void SoPlex2::changeLhsReal(int i, Real lhs)
   {
      assert(_realLP != 0);

      _changeLhsReal(i, lhs);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeLhs(i, lhs);

      _invalidateSolution();
   }



   /// changes right-hand side vector to \p rhs
   void SoPlex2::changeRhsReal(const VectorReal& rhs)
   {
      assert(_realLP != 0);

      _changeRhsReal(rhs);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeRhs(DVectorRational(rhs));

      _invalidateSolution();
   }



   /// changes right-hand side of row \p i to \p rhs
   void SoPlex2::changeRhsReal(int i, Real rhs)
   {
      assert(_realLP != 0);

      _changeRhsReal(i, rhs);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeRhs(i, rhs);

      _invalidateSolution();
   }



   /// changes left- and right-hand side vectors
   void SoPlex2::changeRangeReal(const VectorReal& lhs, const VectorReal& rhs)
   {
      assert(_realLP != 0);

      _changeRangeReal(lhs, rhs);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeRange(DVectorRational(lhs), DVectorRational(rhs));

      _invalidateSolution();
   }



   /// changes left- and right-hand side of row \p i
   void SoPlex2::changeRangeReal(int i, Real lhs, Real rhs)
   {
      assert(_realLP != 0);

      _changeRangeReal(i,lhs, rhs);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeRange(i, lhs, rhs);

      _invalidateSolution();
   }



   /// replaces column \p i with \p lpcol
   void SoPlex2::changeColReal(int i, const LPColReal& lpcol)
   {
      assert(_realLP != 0);

      _changeColReal(i, lpcol);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeCol(i, lpcol);

      _invalidateSolution();
   }



   /// changes vector of lower bounds to \p lower
   void SoPlex2::changeLowerReal(const VectorReal& lower)
   {
      assert(_realLP != 0);

      _changeLowerReal(lower);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeLower(DVectorRational(lower));

      _invalidateSolution();
   }



   /// changes lower bound of column i to \p lower
   void SoPlex2::changeLowerReal(int i, Real lower)
   {
      assert(_realLP != 0);

      _changeLowerReal(i, lower);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeLower(i, lower);

      _invalidateSolution();
   }



   /// changes vector of upper bounds to \p upper
   void SoPlex2::changeUpperReal(const VectorReal& upper)
   {
      assert(_realLP != 0);

      _changeUpperReal(upper);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeUpper(DVectorRational(upper));

      _invalidateSolution();
   }



   /// changes \p i 'th upper bound to \p upper
   void SoPlex2::changeUpperReal(int i, Real upper)
   {
      assert(_realLP != 0);

      _changeUpperReal(i, upper);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeUpper(i, upper);

      _invalidateSolution();
   }



   /// changes vectors of column bounds to \p lower and \p upper
   void SoPlex2::changeBoundsReal(const VectorReal& lower, const VectorReal& upper)
   {
      assert(_realLP != 0);

      _changeBoundsReal(lower, upper);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeBounds(DVectorRational(lower), DVectorRational(upper));

      _invalidateSolution();
   }



   /// changes bounds of column \p i to \p lower and \p upper
   void SoPlex2::changeBoundsReal(int i, Real lower, Real upper)
   {
      assert(_realLP != 0);

      _changeBoundsReal(i, lower, upper);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeBounds(i, lower, upper);

      _invalidateSolution();
   }



   /// changes objective function vector to \p obj
   void SoPlex2::changeObjReal(const VectorReal& obj)
   {
      assert(_realLP != 0);

      _realLP->changeObj(obj);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeObj(DVectorRational(obj));

      _invalidateSolution();
   }



   /// changes objective coefficient of column i to \p obj
   void SoPlex2::changeObjReal(int i, Real obj)
   {
      assert(_realLP != 0);

      _realLP->changeObj(i, obj);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeObj(i, obj);

      _invalidateSolution();
   }



   /// changes matrix entry in row \p i and column \p j to \p val
   void SoPlex2::changeElementReal(int i, int j, Real val)
   {
      assert(_realLP != 0);

      _changeElementReal(i, j, val);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->changeElement(i, j, val);

      _invalidateSolution();
   }



   /// removes row \p i
   void SoPlex2::removeRowReal(int i)
   {
      assert(_realLP != 0);

      _removeRowReal(i);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->removeRow(i);

      _invalidateSolution();
   }



   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsReal()
   void SoPlex2::removeRowsReal(int perm[])
   {
      assert(_realLP != 0);

      _removeRowsReal(perm);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->removeRows(perm);

      _invalidateSolution();
   }



   /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsReal() may be passed
   /// as buffer memory
   void SoPlex2::removeRowsReal(int idx[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numRowsReal());
         _idxToPerm(idx, n, p.get_ptr(), numRowsReal());
         SoPlex2::removeRowsReal(p.get_ptr());
      }
      else
      {
         _idxToPerm(idx, n, perm, numRowsReal());
         SoPlex2::removeRowsReal(perm);
      }
   }



   /// removes rows \p start to \p end including both; an array \p perm of size #numRowsReal() may be passed as buffer
   /// memory
   void SoPlex2::removeRowRangeReal(int start, int end, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numRowsReal());
         _rangeToPerm(start, end, p.get_ptr(), numRowsReal());
         SoPlex2::removeRowsReal(p.get_ptr());
      }
      else
      {
         _rangeToPerm(start, end, perm, numRowsReal());
         SoPlex2::removeRowsReal(perm);
      }
   }



   /// removes column i
   void SoPlex2::removeColReal(int i)
   {
      assert(_realLP != 0);

      _removeColReal(i);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->removeCol(i);

      _invalidateSolution();
   }



   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsReal()
   void SoPlex2::removeColsReal(int perm[])
   {
      assert(_realLP != 0);

      _removeColsReal(perm);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->removeCols(perm);

      _invalidateSolution();
   }



   /// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsReal() may be
   /// passed as buffer memory
   void SoPlex2::removeColsReal(int idx[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numColsReal());
         _idxToPerm(idx, n, p.get_ptr(), numColsReal());
         SoPlex2::removeColsReal(p.get_ptr());
      }
      else
      {
         _idxToPerm(idx, n, perm, numColsReal());
         SoPlex2::removeColsReal(perm);
      }
   }



   /// removes columns \p start to \p end including both; an array \p perm of size #numColsReal() may be passed as
   /// buffer memory
   void SoPlex2::removeColRangeReal(int start, int end, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numColsReal());
         _rangeToPerm(start, end, p.get_ptr(), numColsReal());
         SoPlex2::removeColsReal(p.get_ptr());
      }
      else
      {
         _rangeToPerm(start, end, perm, numColsReal());
         SoPlex2::removeColsReal(perm);
      }
   }



   /// clears the LP
   void SoPlex2::clearLPReal()
   {
      assert(_realLP != 0);

      _realLP->clear();
      _hasBasis = false;

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _rationalLP->clear();

      _invalidateSolution();
   }



   /// synchronizes real LP with rational LP, i.e., copies (rounded) rational LP into real LP, if sync mode is manual
   void SoPlex2::syncLPReal()
   {
      assert(_isConsistent());

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_MANUAL )
         _syncLPReal();
   }



   /// adds a single row
   void SoPlex2::addRowRational(const LPRowRational& lprow)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->addRow(lprow);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _addRowReal(lprow);

      _invalidateSolution();
   }



   /// adds multiple rows
   void SoPlex2::addRowsRational(const LPRowSetRational& lprowset)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->addRows(lprowset);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _addRowsReal(lprowset);

      _invalidateSolution();
   }



   /// adds a single column
   void SoPlex2::addColRational(const LPColRational& lpcol)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->addCol(lpcol);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _addColReal(lpcol);

      _invalidateSolution();
   }



   /// adds multiple columns
   void SoPlex2::addColsRational(const LPColSetRational& lpcolset)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->addCols(lpcolset);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _addColsReal(lpcolset);

      _invalidateSolution();
   }



   /// replaces row \p i with \p lprow
   void SoPlex2::changeRowRational(int i, const LPRowRational& lprow)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeRow(i, lprow);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeRowReal(i, lprow);

      _invalidateSolution();
   }



   /// changes left-hand side vector for constraints to \p lhs
   void SoPlex2::changeLhsRational(const VectorRational& lhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeLhs(lhs);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeLhsReal(DVectorReal(lhs));

      _invalidateSolution();
   }



   /// changes left-hand side of row \p i to \p lhs
   void SoPlex2::changeLhsRational(int i, Rational lhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeLhs(i, lhs);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeLhsReal(i, Real(lhs));

      _invalidateSolution();
   }



   /// changes right-hand side vector to \p rhs
   void SoPlex2::changeRhsRational(const VectorRational& rhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeRhs(rhs);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeRhsReal(DVectorReal(rhs));

      _invalidateSolution();
   }



   /// changes right-hand side of row \p i to \p rhs
   void SoPlex2::changeRhsRational(int i, Rational rhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeRhs(i, rhs);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeRhsReal(i, Real(rhs));

      _invalidateSolution();
   }



   /// changes left- and right-hand side vectors
   void SoPlex2::changeRangeRational(const VectorRational& lhs, const VectorRational& rhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeRange(lhs, rhs);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeRangeReal(DVectorReal(lhs), DVectorReal(rhs));

      _invalidateSolution();
   }



   /// changes left- and right-hand side of row \p i
   void SoPlex2::changeRangeRational(int i, Rational lhs, Rational rhs)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeRange(i, lhs, rhs);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeRangeReal(i, Real(lhs), Real(rhs));

      _invalidateSolution();
   }



   /// replaces column \p i with \p lpcol
   void SoPlex2::changeColRational(int i, const LPColRational& lpcol)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeCol(i, lpcol);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeColReal(i, lpcol);

      _invalidateSolution();
   }



   /// changes vector of lower bounds to \p lower
   void SoPlex2::changeLowerRational(const VectorRational& lower)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeLower(lower);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeLowerReal(DVectorReal(lower));

      _invalidateSolution();
   }



   /// changes lower bound of column i to \p lower
   void SoPlex2::changeLowerRational(int i, Rational lower)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeLower(i, lower);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeLowerReal(i, Real(lower));

      _invalidateSolution();
   }



   /// changes vector of upper bounds to \p upper
   void SoPlex2::changeUpperRational(const VectorRational& upper)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeUpper(upper);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeUpperReal(DVectorReal(upper));

      _invalidateSolution();
   }



   /// changes \p i 'th upper bound to \p upper
   void SoPlex2::changeUpperRational(int i, Rational upper)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeUpper(i, upper);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeUpperReal(i, Real(upper));

      _invalidateSolution();
   }



   /// changes vectors of column bounds to \p lower and \p upper
   void SoPlex2::changeBoundsRational(const VectorRational& lower, const VectorRational& upper)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeBounds(lower, upper);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeBoundsReal(DVectorReal(lower), DVectorReal(upper));

      _invalidateSolution();
   }



   /// changes bounds of column \p i to \p lower and \p upper
   void SoPlex2::changeBoundsRational(int i, Rational lower, Rational upper)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeBounds(i, lower, upper);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeBoundsReal(i, Real(lower), Real(upper));

      _invalidateSolution();
   }



   /// changes objective function vector to \p obj
   void SoPlex2::changeObjRational(const VectorRational& obj)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeObj(obj);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _realLP->changeObj(DVectorReal(obj));

      _invalidateSolution();
   }



   /// changes objective coefficient of column i to \p obj
   void SoPlex2::changeObjRational(int i, Rational obj)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeObj(i, obj);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _realLP->changeObj(i, Real(obj));

      _invalidateSolution();
   }



   /// changes matrix entry in row \p i and column \p j to \p val
   void SoPlex2::changeElementRational(int i, int j, Rational val)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->changeElement(i, j, val);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _changeElementReal(i, j, Real(val));

      _invalidateSolution();
   }



   /// removes row \p i
   void SoPlex2::removeRowRational(int i)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->removeRow(i);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _removeRowReal(i);

      _invalidateSolution();
   }



   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the new
   /// index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsRational()
   void SoPlex2::removeRowsRational(int perm[])
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->removeRows(perm);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _removeRowsReal(perm);

      _invalidateSolution();
   }



   /// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsRational() may be
   /// passed as buffer memory
   void SoPlex2::removeRowsRational(int idx[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numRowsRational());
         _idxToPerm(idx, n, p.get_ptr(), numRowsRational());
         SoPlex2::removeRowsRational(p.get_ptr());
      }
      else
      {
         _idxToPerm(idx, n, perm, numRowsRational());
         SoPlex2::removeRowsRational(perm);
      }
   }



   /// removes rows \p start to \p end including both; an array \p perm of size #numRowsRational() may be passed as
   /// buffer memory
   void SoPlex2::removeRowRangeRational(int start, int end, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numRowsRational());
         _rangeToPerm(start, end, p.get_ptr(), numRowsRational());
         SoPlex2::removeRowsRational(p.get_ptr());
      }
      else
      {
         _rangeToPerm(start, end, perm, numRowsRational());
         SoPlex2::removeRowsRational(perm);
      }
   }



   /// removes column i
   void SoPlex2::removeColRational(int i)
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->removeCol(i);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _removeColReal(i);

      _invalidateSolution();
   }



   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsRational()
   void SoPlex2::removeColsRational(int perm[])
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->removeCols(perm);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
         _removeColsReal(perm);

      _invalidateSolution();
   }



   /// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsRational() may be
   /// passed as buffer memory
   void SoPlex2::removeColsRational(int idx[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numColsRational());
         _idxToPerm(idx, n, p.get_ptr(), numColsRational());
         SoPlex2::removeColsRational(p.get_ptr());
      }
      else
      {
         _idxToPerm(idx, n, perm, numColsRational());
         SoPlex2::removeColsRational(perm);
      }
   }



   /// removes columns \p start to \p end including both; an array \p perm of size #numColsRational() may be passed as
   /// buffer memory
   void SoPlex2::removeColRangeRational(int start, int end, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numColsRational());
         _rangeToPerm(start, end, p.get_ptr(), numColsRational());
         SoPlex2::removeColsRational(p.get_ptr());
      }
      else
      {
         _rangeToPerm(start, end, perm, numColsRational());
         SoPlex2::removeColsRational(perm);
      }
   }



   /// clears the LP
   void SoPlex2::clearLPRational()
   {
      assert(_rationalLP != 0);

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return;

      _rationalLP->clear();

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
      {
         _realLP->clear();
         _hasBasis = false;
      }

      _invalidateSolution();
   }



   /// synchronizes rational LP with real LP, i.e., copies real LP to rational LP, if sync mode is manual
   void SoPlex2::syncLPRational()
   {
      assert(_isConsistent());

      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_MANUAL )
         _syncLPRational();
   }



   /// solves the LP
   SPxSolver::Status SoPlex2::solve()
   {
      assert(_isConsistent());

      // clear statistics
      _statistics->clearSolvingData();

      // the solution is no longer valid
      _invalidateSolution();

      // decide whether to solve the rational LP with iterative refinement or call the standard floating-point solver
      if( Rational::precision() < INT_MAX
          || intParam(SoPlex2::SOLVEMODE) == SOLVEMODE_REAL || (intParam(SoPlex2::SOLVEMODE) == SOLVEMODE_AUTO
            && GE(Real(rationalParam(SoPlex2::FEASTOL)), 1e-9) && GE(Real(rationalParam(SoPlex2::OPTTOL)), 1e-9)) )
      {
         _solveReal();
      }
      else if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
      {
         _syncLPRational();
         _solveRational();
      }
      else if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_MANUAL )
      {
#ifdef ENABLE_ADDITIONAL_CHECKS
         assert(areLPsInSync(true, true, false));
#else
         assert(areLPsInSync(true, false, false));
#endif

         // store current real LP
         SPxLPReal realLP(*_realLP);

         // call rational LP solving with iterative refinement
         _solveRational();

         // restore real LP in order to ensure that we use the same rounding
         assert(_isRealLPLoaded);
         if( _hasBasis )
         {
            assert(_solver.basis().status() > SPxBasis::NO_PROBLEM);
            _basisStatusRows.reSize(_solver.nRows());
            _basisStatusCols.reSize(_solver.nCols());
            _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
         }

         _solver.loadLP(realLP);

         if( _hasBasis )
            _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else
      {
#ifdef ENABLE_ADDITIONAL_CHECKS
         assert(areLPsInSync(true, true, false));
#else
         assert(areLPsInSync(true, false, false));
#endif

         _solveRational();
      }

      return status();
   };



   /// returns the current solver status
   SPxSolver::Status SoPlex2::status() const
   {
      return _status;
   }



   /// is a primal feasible solution available?
   bool SoPlex2::hasPrimal() const
   {
      return (_hasSolReal && _solReal.hasPrimal()) || (_hasSolRational && _solRational.hasPrimal());
   }



   /// is a primal unbounded ray available?
   bool SoPlex2::hasPrimalRay() const
   {
      return (_hasSolReal && _solReal.hasPrimalRay()) || (_hasSolRational && _solRational.hasPrimalRay());
   }



   /// is a dual feasible solution available?
   bool SoPlex2::hasDual() const
   {
      return (_hasSolReal && _solReal.hasDual()) || (_hasSolRational && _solRational.hasDual());
   }



   /// is Farkas proof of infeasibility available?
   bool SoPlex2::hasDualFarkas() const
   {
      return (_hasSolReal && _solReal.hasDualFarkas()) || (_hasSolRational && _solRational.hasDualFarkas());
   }



   /// returns the objective value if a primal solution is available
   ///@todo buffer objective value if computed once
   Real SoPlex2::objValueReal()
   {
      assert(OBJSENSE_MAXIMIZE == 1);
      assert(OBJSENSE_MINIMIZE == -1);

      if( status() == SPxSolver::UNBOUNDED )
         return realParam(SoPlex2::INFTY) * intParam(SoPlex2::OBJSENSE);
      else if( status() == SPxSolver::INFEASIBLE )
         return -realParam(SoPlex2::INFTY) * intParam(SoPlex2::OBJSENSE);
      else if( hasPrimal() )
      {
         _syncRealSolution();
         return (_solReal._primal * maxObjReal()) * intParam(SoPlex2::OBJSENSE);
      }
      else
         return 0;
   }



   /// gets the primal solution vector if available; returns true on success
   bool SoPlex2::getPrimalReal(VectorReal& vector)
   {
      if( hasPrimal() && vector.dim() >= numColsReal() )
      {
         _syncRealSolution();
         _solReal.getPrimal(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the vector of slack values if available; returns true on success
   bool SoPlex2::getSlacksReal(VectorReal& vector)
   {
      if( hasPrimal() && vector.dim() >= numRowsReal() )
      {
         _syncRealSolution();
         _solReal.getSlacks(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the primal ray if available; returns true on success
   bool SoPlex2::getPrimalRayReal(VectorReal& vector)
   {
      if( hasPrimalRay() && vector.dim() >= numColsReal() )
      {
         _syncRealSolution();
         _solReal.getPrimalRay(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the dual solution vector if available; returns true on success
   bool SoPlex2::getDualReal(VectorReal& vector)
   {
      if( hasDual() && vector.dim() >= numRowsReal() )
      {
         _syncRealSolution();
         _solReal.getDual(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the vector of reduced cost values if available; returns true on success
   bool SoPlex2::getRedCostReal(VectorReal& vector)
   {
      if( hasDual() && vector.dim() >= numColsReal() )
      {
         _syncRealSolution();
         _solReal.getRedCost(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the Farkas proof if available; returns true on success
   bool SoPlex2::getDualFarkasReal(VectorReal& vector)
   {
      if( hasDualFarkas() && vector.dim() >= numRowsReal() )
      {
         _syncRealSolution();
         _solReal.getDualFarkas(vector);
         return true;
      }
      else
         return false;
   }



   /// gets violation of bounds by given primal solution
   void SoPlex2::getBoundViolationReal(VectorReal& primal, Real& maxviol, Real& sumviol) const
   {
      assert(primal.dim() >= numColsReal());

      if( primal.dim() < numColsReal() )
      {
         maxviol = realParam(SoPlex2::INFTY);
         sumviol = realParam(SoPlex2::INFTY);
         return;
      }

      maxviol = 0.0;
      sumviol = 0.0;

      for( int i = numColsReal() - 1; i >= 0; i-- )
      {
         Real viol = lowerReal(i) - primal[i];
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }

         viol = primal[i] - upperReal(i);
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }
      }
   }



   /// gets violation of constraints
   void SoPlex2::getConstraintViolationReal(VectorReal& primal, Real& maxviol, Real& sumviol) const
   {
      assert(primal.dim() >= numColsReal());

      if( primal.dim() < numColsReal() )
      {
         maxviol = realParam(SoPlex2::INFTY);
         sumviol = realParam(SoPlex2::INFTY);
         return;
      }

      DVectorReal activity = _realLP->computePrimalActivity(primal);
      maxviol = 0.0;
      sumviol = 0.0;

      for( int i = numRowsReal() - 1; i >= 0; i-- )
      {
         Real viol = lhsReal(i) - activity[i];
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }

         viol = activity[i] - rhsReal(i);
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }
      }
   }



   /// gets violation of slacks
   ///@todo implement
   void SoPlex2::getSlackViolationReal(Real& maxviol, Real& sumviol) const
   {
      maxviol = 0.0;
      sumviol = 0.0;
   }



   /// gets violation of reduced costs
   ///@todo implement
   void SoPlex2::getRedCostViolationReal(Real& maxviol, Real& sumviol) const
   {
      maxviol = 0.0;
      sumviol = 0.0;
   }



   /// returns the objective value if a primal solution is available
   Rational SoPlex2::objValueRational()
   {
      assert(OBJSENSE_MAXIMIZE == 1);
      assert(OBJSENSE_MINIMIZE == -1);

      if( status() == SPxSolver::UNBOUNDED )
         return Rational(realParam(SoPlex2::INFTY) * intParam(SoPlex2::OBJSENSE));
      else if( status() == SPxSolver::INFEASIBLE )
         return Rational(-realParam(SoPlex2::INFTY) * intParam(SoPlex2::OBJSENSE));
      else if( hasPrimal() )
      {
         _syncRationalSolution();
         return (_solRational._primal * maxObjRational()) * intParam(SoPlex2::OBJSENSE);
      }
      else
         return 0;
   }



   /// gets the primal solution vector if available; returns true on success
   bool SoPlex2::getPrimalRational(VectorRational& vector)
   {
      if( hasPrimal() && vector.dim() >= numColsRational() )
      {
         _syncRationalSolution();
         _solRational.getPrimal(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the vector of slack values if available; returns true on success
   bool SoPlex2::getSlacksRational(VectorRational& vector)
   {
      if( hasPrimal() && vector.dim() >= numRowsRational() )
      {
         _syncRationalSolution();
         _solRational.getSlacks(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the primal ray if LP is unbounded; returns true on success
   bool SoPlex2::getPrimalRayRational(VectorRational& vector)
   {
      if( hasPrimalRay() && vector.dim() >= numColsRational() )
      {
         _syncRationalSolution();
         _solRational.getPrimalRay(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the dual solution vector if available; returns true on success
   bool SoPlex2::getDualRational(VectorRational& vector)
   {
      if( hasDual() && vector.dim() >= numRowsRational() )
      {
         _syncRationalSolution();
         _solRational.getDual(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the vector of reduced cost values if available; returns true on success
   bool SoPlex2::getRedCostRational(VectorRational& vector)
   {
      if( hasDual() && vector.dim() >= numColsRational() )
      {
         _syncRationalSolution();
         _solRational.getRedCost(vector);
         return true;
      }
      else
         return false;
   }



   /// gets the Farkas proof if LP is infeasible; returns true on success
   bool SoPlex2::getDualFarkasRational(VectorRational& vector)
   {
      if( hasDualFarkas() && vector.dim() >= numRowsRational() )
      {
         _syncRationalSolution();
         _solRational.getDualFarkas(vector);
         return true;
      }
      else
         return false;
   }



   /// gets violation of bounds by given primal solution
   void SoPlex2::getBoundViolationRational(VectorRational& primal, Rational& maxviol, Rational& sumviol) const
   {
      assert(primal.dim() >= numColsRational());

      if( primal.dim() < numColsRational() )
      {
         maxviol = Rational(realParam(SoPlex2::INFTY));
         sumviol = Rational(realParam(SoPlex2::INFTY));
         return;
      }

      maxviol = 0;
      sumviol = 0;

      for( int i = numColsRational() - 1; i >= 0; i-- )
      {
         Rational viol = lowerRational(i) - primal[i];
         if( viol > 0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }

         viol = primal[i] - upperRational(i);
         if( viol > 0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }
      }
   }



   /// gets violation of constraints
   void SoPlex2::getConstraintViolationRational(VectorRational& primal, Rational& maxviol, Rational& sumviol) const
   {
      assert(primal.dim() >= numColsRational());

      if( primal.dim() < numColsRational() )
      {
         maxviol = Rational(realParam(SoPlex2::INFTY));
         sumviol = Rational(realParam(SoPlex2::INFTY));
         return;
      }

      DVectorRational activity = _rationalLP->computePrimalActivity(primal);
      maxviol = 0;
      sumviol = 0;

      for( int i = numRowsRational() - 1; i >= 0; i-- )
      {
         Rational viol = lhsRational(i) - activity[i];
         if( viol > 0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }

         viol = activity[i] - rhsRational(i);
         if( viol > 0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }
      }
   }



   /// gets violation of slacks
   ///@todo implement
   void SoPlex2::getSlackViolationRational(Rational& maxviol, Rational& sumviol) const
   {
      maxviol = 0;
      sumviol = 0;
   }



   /// gets violation of reduced costs
   ///@todo implement
   void SoPlex2::getRedCostViolationRational(Rational& maxviol, Rational& sumviol) const
   {
      maxviol = 0;
      sumviol = 0;
   }



   /// is an advanced starting basis available?
   bool SoPlex2::hasBasis() const
   {
      return _hasBasis;
   }



   /// returns the current basis status
   SPxBasis::SPxStatus SoPlex2::basisStatus() const
   {
      if( !hasBasis() )
         return SPxBasis::NO_PROBLEM;
      else if( status() == SPxSolver::OPTIMAL )
         return SPxBasis::OPTIMAL;
      else if( status() == SPxSolver::UNBOUNDED )
         return SPxBasis::UNBOUNDED;
      else if( status() == SPxSolver::INFEASIBLE )
         return SPxBasis::INFEASIBLE;
      else if( hasPrimal() )
         return SPxBasis::PRIMAL;
      else if( hasDual() )
         return SPxBasis::DUAL;
      else
         return SPxBasis::REGULAR;
   }



   /// returns basis status for a single row
   SPxSolver::VarStatus SoPlex2::basisRowStatus(int row) const
   {
      assert(row >= 0);
      assert(row < numRowsReal());

      // if no basis is available, return slack basis; if index is out of range, return basic status as for a newly
      // added row
      if( !hasBasis() || row < 0 || row >= numRowsReal() )
         return SPxSolver::BASIC;
      // if the real LP is loaded, ask solver
      else if( _isRealLPLoaded )
      {
         return _solver.getBasisRowStatus(row);
      }
      // if the real LP is not loaded, the basis is stored in the basis arrays of this class
      else
      {
         assert(row < _basisStatusRows.size());
         return _basisStatusRows[row];
      }
   }



   /// returns basis status for a single column
   SPxSolver::VarStatus SoPlex2::basisColStatus(int col) const
   {
      assert(col >= 0);
      assert(col < numColsReal());

      // if index is out of range, return nonbasic status as for a newly added unbounded column
      if( col < 0 || col >= numColsReal() )
      {
         return SPxSolver::ZERO;
      }
      // if no basis is available, return slack basis
      else if( !hasBasis() )
      {
         if( lowerReal(col) > -realParam(SoPlex2::INFTY) )
            return SPxSolver::ON_LOWER;
         else if( upperReal(col) < realParam(SoPlex2::INFTY) )
            return SPxSolver::ON_UPPER;
         else
            return SPxSolver::ZERO;
      }
      // if the real LP is loaded, ask solver
      else if( _isRealLPLoaded )
      {
         return _solver.getBasisColStatus(col);
      }
      // if the real LP is not loaded, the basis is stored in the basis arrays of this class
      else
      {
         assert(col < _basisStatusCols.size());
         return _basisStatusCols[col];
      }
   }



   /// gets current basis
   void SoPlex2::getBasis(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[]) const
   {
      // if no basis is available, return slack basis
      if( !hasBasis() )
      {
         for( int i = numRowsReal() - 1; i >= 0; i-- )
            rows[i] = SPxSolver::BASIC;

         for( int i = numColsReal() - 1; i >= 0; i-- )
         {
            if( lowerReal(i) > -realParam(SoPlex2::INFTY) )
               cols[i] = SPxSolver::ON_LOWER;
            else if( upperReal(i) < realParam(SoPlex2::INFTY) )
               cols[i] = SPxSolver::ON_UPPER;
            else
               cols[i] = SPxSolver::ZERO;
         }
      }
      // if the real LP is loaded, ask solver
      else if( _isRealLPLoaded )
      {
         (void)_solver.getBasis(rows, cols);
      }
      // if the real LP is not loaded, the basis is stored in the basis arrays of this class
      else
      {
         assert(numRowsReal() == _basisStatusRows.size());
         assert(numColsReal() == _basisStatusCols.size());

         for( int i = numRowsReal() - 1; i >= 0; i-- )
            rows[i] = _basisStatusRows[i];

         for( int i = numColsReal() - 1; i >= 0; i-- )
            cols[i] = _basisStatusCols[i];
      }
   }



   /// returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m
   void SoPlex2::getBasisInd(int* bind) const
   {
      // if no basis is available, return slack basis
      if( !hasBasis() )
      {
         for( int i = 0; i < numRowsReal(); ++i )
            bind[i] = -1 - i;
      }
      // if the real LP is not loaded, the basis is stored in the basis arrays of this class
      else if( !_isRealLPLoaded )
      {
         int k = 0;

         assert(numRowsReal() == _basisStatusRows.size());
         assert(numColsReal() == _basisStatusCols.size());

         for( int i = 0; i < numRowsReal(); ++i )
         {
            if( _basisStatusRows[i] == SPxSolver::BASIC )
            {
               bind[k] = -1 - i;
               k++;
            }
         }

         for( int j = 0; j < numColsReal(); ++j )
         {
            if( _basisStatusCols[j] == SPxSolver::BASIC )
            {
               bind[k] = j;
               k++;
            }
         }

         assert(k == numRowsReal());
      }
      // if the real LP is loaded, the basis is stored in the solver and we need to distinguish between column and row
      // representation; ask the solver itself which representation it has, since the REPRESENTATION parameter of this
      // class might be set to automatic
      else if( _solver.rep() == SPxSolver::COLUMN )
      {
         for( int i = 0; i < numRowsReal(); ++i )
         {
            SPxId id = _solver.basis().baseId(i);
            bind[i] = (id.isSPxColId() ? _solver.number(id) : - 1 - _solver.number(id));
         }
      }
      // for row representation, return the complement of the basis; for this, we need to loop through all rows and columns
      else
      {
         assert(_solver.rep() == SPxSolver::ROW);

         int k = 0;

         for( int i = 0; i < numRowsReal(); ++i )
         {
            if( !_solver.isRowBasic(i) )
            {
               bind[k] = -1 - i;
               k++;
            }
         }

         for( int j = 0; j < numColsReal(); ++j )
         {
            if( !_solver.isColBasic(j) )
            {
               bind[k] = j;
               k++;
            }
         }

         assert(k == numRowsReal());
      }
   }



   /// computes an estimated condition number for the current basis matrix using the power method; returns true on success
   bool SoPlex2::getEstimatedCondition(Real& condition)
   {
      _ensureRealLPLoaded();
      if( !_isRealLPLoaded )
         return false;

      if( _solver.basis().status() == SPxBasis::NO_PROBLEM )
         return false;

      condition = _solver.basis().getEstimatedCondition();

      return true;
   }

   /// computes the exact condition number for the current basis matrix using the power method; returns true on success
   bool SoPlex2::getExactCondition(Real& condition)
   {
      _ensureRealLPLoaded();
      if( !_isRealLPLoaded )
         return false;

      if( _solver.basis().status() == SPxBasis::NO_PROBLEM )
         return false;

      condition = _solver.basis().getExactCondition();

      return true;
   }

   /// computes row r of basis inverse; returns true on success
   ///@todo use VectorReal for coef
   bool SoPlex2::getBasisInverseRowReal(int r, Real* coef)
   {
      assert(r >= 0);
      assert(r < numRowsReal());
      assert(coef != 0);

      if( !hasBasis() || r < 0 || r >= numRowsReal() )
         return false;

      _ensureRealLPLoaded();

      if( !_isRealLPLoaded )
         return false;

      // we need to distinguish between column and row representation; ask the solver itself which representation it
      // has, since the REPRESENTATION parameter of this class might be set to automatic
      if( _solver.rep() == SPxSolver::COLUMN )
      {
         SSVectorReal x(numRowsReal());
         try
         {
            _solver.basis().coSolve(x, _solver.unitVector(r));
         }
         catch( SPxException E )
         {
            MSG_ERROR( spxout << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
         }
         // copy sparse data to dense result vector based on coef array
         VectorReal y(numRowsReal(), coef);
         y = x;
      }
      else
      {
         assert(_solver.rep() == SPxSolver::ROW);

         // @todo should rhs be a reference?
         DSVector rhs(numColsReal());
         SSVector y(numColsReal());
         int* bind = 0;
         int index;

         // get ordering of column basis matrix
         spx_alloc(bind, numRowsReal());
         getBasisInd(bind);

         // get vector corresponding to requested index r
         index = bind[r];

         // r corresponds to a row vector
         if( index < 0 )
         {
            // transform index to actual row index
            index = -index - 1;

            // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index >= 0);
            assert(index < numRowsReal());
            assert(!_solver.isRowBasic(index));

            // get row vector
            rhs = _solver.rowVector(index);
            rhs *= -1.0;
         }
         // r corresponds to a column vector
         else
         {
            // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index < numColsReal());
            assert(!_solver.isColBasic(index));

            // get unit vector
            rhs = UnitVector(index);
         }

         // solve system "y B = rhs", where B is the row basis matrix
         try
         {
            _solver.basis().solve(y, rhs);
         }
         catch( SPxException E )
         {
            MSG_ERROR( spxout << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
         }

         // initialize result vector x as zero
         memset(coef, 0, numRowsReal() * sizeof(Real));

         // add nonzero entries
         for( int i = 0; i < numColsReal(); ++i )
         {
            SPxId id = _solver.basis().baseId(i);

            if( id.isSPxRowId() )
            {
               assert(_solver.number(id) >= 0);
               assert(_solver.number(id) < numRowsReal());
               assert(bind[r] >= 0 || _solver.number(id) != index);

               coef[_solver.number(id)] = y[i];
            }
         }

         // if r corresponds to a row vector, we have to add a 1 at position r
         if( bind[r] < 0 )
         {
            assert(coef[index] == 0.0);
            coef[index] = 1.0;
         }

         // free memory
         spx_free(bind);
      }

      return true;
   }



   /// computes column c of basis inverse; returns true on success
   ///@todo use VectorReal for coef
   bool SoPlex2::getBasisInverseColReal(int c, Real* coef)
   {
      assert(c >= 0);
      assert(c < numRowsReal());
      assert(coef != 0);

      if( !hasBasis() || c < 0 || c >= numRowsReal() )
         return false;

      _ensureRealLPLoaded();

      if( !_isRealLPLoaded )
         return false;

      // we need to distinguish between column and row representation; ask the solver itself which representation it
      // has, since the REPRESENTATION parameter of this class might be set to automatic
      if( _solver.rep() == SPxSolver::COLUMN )
      {
         SSVectorReal x(numColsReal());
         try
         {
            _solver.basis().solve(x, _solver.unitVector(c));
         }
         catch( SPxException E )
         {
            MSG_ERROR( spxout << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
         }
         // copy sparse data to dense result vector based on coef array
         VectorReal y(numColsReal(), coef);
         y = x;
      }
      else
      {
         assert(_solver.rep() == SPxSolver::ROW);

         // @todo should rhs be a reference?
         DSVectorReal rhs(numColsReal());
         SSVectorReal y(numColsReal());
         int* bind = 0;
         int index;

         // get ordering of column basis matrix
         spx_alloc(bind, numRowsReal());
         getBasisInd(bind);

         // get vector corresponding to requested index r
         index = bind[c];

         // c corresponds to a row vector
         if( index < 0 )
         {
            // transform index to actual row index
            index = -index - 1;

            // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index >= 0);
            assert(index < numRowsReal());
            assert(!_solver.isRowBasic(index));

            // get row vector
            rhs = _solver.rowVector(index);
            rhs *= -1.0;
         }
         // c corresponds to a column vector
         else
         {
            // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index < numColsReal());
            assert(!_solver.isColBasic(index));

            // get unit vector
            rhs = UnitVector(index);
         }

         // solve system "y B = rhs", where B is the row basis matrix
         try
         {
            _solver.basis().coSolve(y, rhs);
         }
         catch( SPxException E )
         {
            MSG_ERROR( spxout << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
         }

         // initialize result vector x as zero
         memset(coef, 0, numRowsReal() * sizeof(Real));

         // add nonzero entries
         for( int i = 0; i < numColsReal(); ++i )
         {
            SPxId id = _solver.basis().baseId(i);

            if( id.isSPxRowId() )
            {
               assert(_solver.number(id) >= 0);
               assert(_solver.number(id) < numRowsReal());
               assert(bind[c] >= 0 || _solver.number(id) != index);

               coef[_solver.number(id)] = y[i];
            }
         }

         // if r corresponds to a row vector, we have to add a 1 at position r
         if( bind[c] < 0 )
         {
            assert(coef[index] == 0.0);
            coef[index] = 1.0;
         }

         // free memory
         spx_free(bind);
      }

      return true;
   }



   /// computes dense solution of basis matrix B * sol = rhs; returns true on success
   ///@todo use VectorReal for rhs and sol
   bool SoPlex2::getBasisInverseTimesVecReal(Real* rhs, Real* sol)
   {
      VectorReal v(numRowsReal(), rhs);
      VectorReal x(numRowsReal(), sol);

      if( !hasBasis() )
         return false;

      _ensureRealLPLoaded();

      if( !_isRealLPLoaded )
         return false;

      // we need to distinguish between column and row representation; ask the solver itself which representation it
      // has, since the REPRESENTATION parameter of this class might be set to automatic; in the column case we can use
      // the existing factorization
      if( _solver.rep() == SPxSolver::COLUMN )
      {
         // solve system "x = B^-1 * A_c" to get c'th column of B^-1 * A
         try
         {
            _solver.basis().solve(x, v);
         }
         catch( SPxException E )
         {
            MSG_ERROR( spxout << "Caught exception <" << E.what() << "> while solving with basis matrix.\n" );
            return false;
         }
      }
      else
      {
         assert(_solver.rep() == SPxSolver::ROW);

         DSVectorReal rowrhs(numColsReal());
         SSVectorReal y(numColsReal());
         int* bind;

         // get ordering of column basis matrix
         spx_alloc(bind, numRowsReal());
         getBasisInd(bind);

         // fill right-hand side for row-based system
         for( int i = 0; i < numColsReal(); ++i )
         {
            SPxId id = _solver.basis().baseId(i);

            if( id.isSPxRowId() )
            {
               assert(_solver.number(id) >= 0);
               assert(_solver.number(id) < numRowsReal());

               rowrhs.add(i, v[_solver.number(id)]);
            }
            else
            {
               assert(rowrhs[i] == 0.0);
            }
         }

         // solve system "B y = rowrhs", where B is the row basis matrix
         try
         {
            _solver.basis().coSolve(y, rowrhs);
         }
         catch( SPxException E )
         {
            MSG_ERROR( spxout << "Caught exception <" << E.what() << "> while solving with basis matrix.\n" );
            return false;
         }

         // fill result w.r.t. order given by bind
         for( int i = 0; i < numRowsReal(); ++i )
         {
            int index;

            index = bind[i];

            if( index < 0 )
            {
               index = -index-1;

               // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
               assert(index >= 0);
               assert(index < numRowsReal());
               assert(!_solver.isRowBasic(index));

               x[i] = v[index] - (rowVectorReal(index) * Vector(numColsReal(), y.get_ptr()));
            }
            else
            {
               // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
               assert(index >= 0);
               assert(index < numColsReal());
               assert(!_solver.isColBasic(index));

               x[i] = y[index];
            }
         }

         // free memory
         spx_free(bind);
      }
      return true;
   }



   /// sets starting basis via arrays of statuses
   void SoPlex2::setBasis(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[])
   {
      if( _isRealLPLoaded )
      {
         assert(numRowsReal() == _solver.nRows());
         assert(numColsReal() == _solver.nCols());

         _solver.setBasis(rows, cols);
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else
      {
         _basisStatusRows.reSize(numRowsReal());
         _basisStatusCols.reSize(numColsReal());

         for( int i = numRowsReal() - 1; i >= 0; i-- )
            _basisStatusRows[i] = rows[i];

         for( int j = numColsReal() - 1; j >= 0; j-- )
            _basisStatusCols[j] = cols[j];

         _hasBasis = true;
      }
   }



   /// clears starting basis
   void SoPlex2::clearBasis()
   {
      if( _isRealLPLoaded )
         _solver.reLoad();

      _hasBasis = false;
   }



   /// number of iterations since last call to solve
   int SoPlex2::numIterations() const
   {
      return _statistics->iterations;
   }



   /// time spent in last call to solve
   Real SoPlex2::solveTime() const
   {
       return _statistics->solvingTime.userTime();
   }



   /// statistical information in form of a string
   std::string SoPlex2::statisticString() const
   {
      std::stringstream s;
      s  << "Factorizations     : " << std::setw(10) << _statistics->luFactorizations << std::endl
         << "  Time spent       : " << std::setw(10) << std::fixed << std::setprecision(2) << _statistics->luFactorizationTime << std::endl
         << "Solves             : " << std::setw(10) << _statistics->luSolves << std::endl
         << "  Time spent       : " << std::setw(10) << _statistics->luSolveTime << std::endl
         << "Solution time      : " << std::setw(10) << std::fixed << std::setprecision(2) << solveTime() << std::endl
         << "Iterations         : " << std::setw(10) << numIterations() << std::endl;

      return s.str();
   }



   /// name of starter
   const char* SoPlex2::getStarterName()
   {
      if( _starter )
         return _starter->getName();
      else
         return "none";
   }



   /// name of simplifier
   const char* SoPlex2::getSimplifierName()
   {
      if( _simplifier )
         return _simplifier->getName();
      else
         return "none";
   }



   /// name of scaling method after simplifier
   const char* SoPlex2::getScalerName()
   {
      if( _scaler )
         return _scaler->getName();
      else
         return "none";
   }



   /// name of currently loaded pricer
   const char* SoPlex2::getPricerName()
   {
      return _solver.pricer()->getName();
   }



   /// name of currently loaded ratiotester
   const char* SoPlex2::getRatiotesterName()
   {
      return _solver.ratiotester()->getName();
   }



   /// reads LP file in LP or MPS format according to READMODE parameter; gets row names, column names, and
   /// integer variables if desired; returns true on success
   bool SoPlex2::readFile(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
   {
      if( intParam(SoPlex2::READMODE) == READMODE_REAL )
         return _readFileReal(filename, rowNames, colNames, intVars);
      else
         return _readFileRational(filename, rowNames, colNames, intVars);
   }

   /// writes real LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer; returns true on success
   bool SoPlex2::writeFileReal(const char* filename, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars) const
   {
      ///@todo implement return value
      _realLP->writeFile(filename, rowNames, colNames, intVars);
      return true;
   }



   /// writes rational LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer; returns true on success
   bool SoPlex2::writeFileRational(const char* filename, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars) const
   {
      if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         return false;
      else
      {
         assert(_rationalLP != 0);
         _rationalLP->writeFile(filename, rowNames, colNames, intVars);

         ///@todo implement return value
         return true;
      }
   }



   /// reads basis information from \p filename and returns true on success; if \p rowNames and \p colNames are \c NULL,
   /// default names are assumed; returns true on success
   bool SoPlex2::readBasisFile(const char* filename, const NameSet* rowNames, const NameSet* colNames)
   {
#if 1
      assert(filename != 0);
      assert(_realLP != 0);

      // start timing
      _statistics->readingTime.start();

      // read
      if( !_isRealLPLoaded )
      {
         assert(_realLP != &_solver);

         _solver.loadLP(*_realLP);
         _realLP->~SPxLPReal();
         spx_free(_realLP);
         _realLP = &_solver;
         _isRealLPLoaded = true;
      }
      _hasBasis = _solver.readBasisFile(filename, rowNames, colNames);
      assert(_hasBasis == (_solver.basis().status() > SPxBasis::NO_PROBLEM));

      // stop timing
      _statistics->readingTime.stop();

      return _hasBasis;
#else
      // this is alternative code for reading bases without the SPxSolver class
      assert(filename != 0);

      // start timing
      _statistics->readingTime.start();

      // read
      spxifstream file(filename);

      if( !file )
         return false;

      // get problem size
      int numRows = numRowsReal();
      int numCols = numColsReal();

      // prepare column names
      const NameSet* colNamesPtr = colNames;
      NameSet* tmpColNames = 0;
      if( colNames == 0 )
      {
         std::stringstream name;

         spx_alloc(tmpColNames);
         tmpColNames = new (tmpColNames) NameSet();
         tmpColNames->reMax(numCols);

         for( int j = 0; j < numCols; ++j )
         {
            name << "x" << j;
            tmpColNames->add(name.str().c_str());
         }

         colNamesPtr = tmpColNames;
      }

      // prepare row names
      const NameSet* rowNamesPtr = rowNames;
      NameSet* tmpRowNames = 0;
      if( rowNamesPtr == 0 )
      {
         std::stringstream name;

         spx_alloc(tmpRowNames);
         tmpRowNames = new (tmpRowNames) NameSet();
         tmpRowNames->reMax(numRows);

         for( int i = 0; i < numRows; ++i )
         {
            name << "C" << i;
            tmpRowNames->add(name.str().c_str());
         }

         rowNamesPtr = tmpRowNames;
      }

      // initialize with default slack basis
      _basisStatusRows.reSize(numRows);
      _basisStatusCols.reSize(numCols);

      for( int i = 0; i < numRows; i++ )
         _basisStatusRows[i] = SPxSolver::BASIC;

      for( int i = 0; i < numCols; i++ )
      {
         if( lowerReal(i) == upperReal(i) )
            _basisStatusCols[i] = SPxSolver::FIXED;
         else if( lowerReal(i) <= double(-realParam(SoPlex2::INFTY)) && upperReal(i) >= double(realParam(SoPlex2::INFTY)) )
            _basisStatusCols[i] = SPxSolver::ZERO;
         else if( lowerReal(i) <= double(-realParam(SoPlex2::INFTY)) )
            _basisStatusCols[i] = SPxSolver::ON_UPPER;
         else
            _basisStatusCols[i] = SPxSolver::ON_LOWER;
      }

      // read basis
      MPSInput mps(file);
      if( mps.readLine() && (mps.field0() != 0) && !strcmp(mps.field0(), "NAME") )
      {
         while( mps.readLine() )
         {
            int c = -1;
            int r = -1;

            if( mps.field0() != 0 && !strcmp(mps.field0(), "ENDATA") )
            {
               mps.setSection(MPSInput::ENDATA);
               break;
            }

            if( mps.field1() == 0 || mps.field2() == 0 )
               break;

            if( (c = colNamesPtr->number(mps.field2())) < 0 )
               break;

            if( *mps.field1() == 'X' )
            {
               if( mps.field3() == 0 || (r = rowNamesPtr->number(mps.field3())) < 0 )
                  break;
            }

            if( !strcmp(mps.field1(), "XU") )
            {
               _basisStatusCols[c] = SPxSolver::BASIC;
               _basisStatusRows[r] = (lhsReal(r) == rhsReal(r))
                  ? SPxSolver::FIXED
                  : SPxSolver::ON_UPPER;
            }
            else if( !strcmp(mps.field1(), "XL") )
            {
               _basisStatusCols[c] = SPxSolver::BASIC;
               _basisStatusRows[r] = (lhsReal(r) == rhsReal(r))
                  ? SPxSolver::FIXED
                  : SPxSolver::ON_LOWER;
            }
            else if( !strcmp(mps.field1(), "UL") )
            {
               _basisStatusCols[c] = SPxSolver::ON_UPPER;
            }
            else if( !strcmp(mps.field1(), "LL") )
            {
               _basisStatusCols[c] = SPxSolver::ON_LOWER;
            }
            else
            {
               mps.syntaxError();
               break;
            }
         }
      }

      if( rowNames == 0 )
      {
         tmpRowNames->~NameSet();
         spx_free(tmpRowNames);
      }

      if( colNames == 0 )
      {
         tmpColNames->~NameSet();
         spx_free(tmpColNames);
      }

      _hasBasis = !mps.hasError();

      // stop timing
      _statistics->readingTime.stop();

      return _hasBasis;
#endif
   }



   /// writes basis information to \p filename; if \p rowNames and \p colNames are \c NULL, default names are used;
   /// returns true on success
   bool SoPlex2::writeBasisFile(const char* filename, const NameSet* rowNames, const NameSet* colNames) const
   {
      assert(filename != 0);

      if( _isRealLPLoaded )
         return _solver.writeBasisFile(filename, rowNames, colNames);
      else
      {
         std::ofstream file(filename);
         if( file == 0 )
            return false;

         file.setf(std::ios::left);
         file << "NAME  " << filename << "\n";

         // do not write basis if there is none
         if( !_hasBasis )
         {
            file << "ENDATA\n";
            return true;
         }

         // start writing
         int numRows = _basisStatusRows.size();
         int numCols = _basisStatusCols.size();
         int row = 0;

         for( int col = 0; col < numCols; col++ )
         {
            assert(_basisStatusCols[col] != SPxSolver::UNDEFINED);

            if( _basisStatusCols[col] == SPxSolver::BASIC )
            {
               // find nonbasic row
               for( ; row < numRows; row++ )
               {
                  assert(_basisStatusRows[row] != SPxSolver::UNDEFINED);
                  if( _basisStatusRows[row] != SPxSolver::BASIC )
                     break;
               }

               assert(row != numRows);

               file << (_basisStatusRows[row] == SPxSolver::ON_UPPER ? " XU " : " XL ");

               file << std::setw(8);
               if( colNames != 0 && colNames->has(col) )
                  file << (*colNames)[col];
               else
                  file << "x" << col;

               file << "       ";
               if( rowNames != 0 && rowNames->has(row) )
                  file << (*rowNames)[row];
               else
                  file << "C" << row;

               file << "\n";
               row++;
            }
            else
            {
               if( _basisStatusCols[col] == SPxSolver::ON_UPPER )
               {
                  file << " UL ";

                  file << std::setw(8);
                  if( colNames != 0 && colNames->has(col) )
                     file << (*colNames)[col];
                  else
                     file << "x" << col;

                  file << "\n";
               }
            }
         }

         file << "ENDATA\n";

#ifndef NDEBUG
         // check that the remaining rows are basic
         for( ; row < numRows; row++ )
         {
            assert(_basisStatusRows[row] == SPxSolver::BASIC);
         }
#endif

         return true;
      }
   }



   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void SoPlex2::writeStateReal(const char* filename, const NameSet* rowNames, const NameSet* colNames) const
   {
      std::string ofname;

      // write parameter settings
      ofname = std::string(filename) + ".set";
      saveSettingsFile(ofname.c_str());

      // write problem in MPS format
      ofname = std::string(filename) + ".mps";
      writeFileReal(ofname.c_str(), rowNames, colNames, 0);

      // write basis
      ofname = std::string(filename) + ".bas";
      writeBasisFile(ofname.c_str(), rowNames, colNames);
   }



   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void SoPlex2::writeStateRational(const char* filename, const NameSet* rowNames, const NameSet* colNames) const
   {
      std::string ofname;

      // write parameter settings
      ofname = std::string(filename) + ".set";
      saveSettingsFile(ofname.c_str());

      // write problem in MPS format
      ofname = std::string(filename) + ".mps";
      writeFileRational(ofname.c_str(), rowNames, colNames, 0);

      // write basis
      ofname = std::string(filename) + ".bas";
      writeBasisFile(ofname.c_str(), rowNames, colNames);
   }



   /// returns boolean parameter value
   bool SoPlex2::boolParam(const BoolParam param) const
   {
      assert(param >= 0);
      assert(param < SoPlex2::BOOLPARAM_COUNT);
      return _currentSettings->_boolParamValues[param];
   }



   /// returns integer parameter value
   int SoPlex2::intParam(const IntParam param) const
   {
      assert(param >= 0);
      assert(param < INTPARAM_COUNT);
      return _currentSettings->_intParamValues[param];
   }



   /// returns real parameter value
   Real SoPlex2::realParam(const RealParam param) const
   {
      assert(param >= 0);
      assert(param < REALPARAM_COUNT);
      return _currentSettings->_realParamValues[param];
   }



   /// returns rational parameter value
   Rational SoPlex2::rationalParam(const RationalParam param) const
   {
      assert(param >= 0);
      assert(param < RATIONALPARAM_COUNT);
      return _currentSettings->_rationalParamValues[param];
   }



   /// returns current parameter settings
   const SoPlex2::Settings& SoPlex2::settings() const
   {
      return *_currentSettings;
   }



   /// sets boolean parameter value; returns true on success
   bool SoPlex2::setBoolParam(const BoolParam param, const bool value, const bool quiet, const bool init)
   {
      assert(param >= 0);
      assert(param < SoPlex2::BOOLPARAM_COUNT);
      assert(init || _isConsistent());

      if( !init && value == boolParam(param) )
         return true;

      switch( param )
      {
      case PARTIAL_PRICING:
         ///@todo activate partial pricing
         break;
      case LIFTING:
         break;
      default:
         return false;
      }

      _currentSettings->_boolParamValues[param] = value;
      return true;
   }



   /// sets integer parameter value; returns true on success
   bool SoPlex2::setIntParam(const IntParam param, const int value, const bool quiet, const bool init)
   {
      assert(param >= 0);
      assert(param < INTPARAM_COUNT);
      assert(init || _isConsistent());

      if( !init && value == intParam(param) )
         return true;

      switch( param )
      {
      // objective sense
      case SoPlex2::OBJSENSE:
         if( value != SoPlex2::OBJSENSE_MAXIMIZE && value != SoPlex2::OBJSENSE_MINIMIZE )
            return false;
         _realLP->changeSense(value == SoPlex2::OBJSENSE_MAXIMIZE ? SPxLPReal::MAXIMIZE : SPxLPReal::MINIMIZE);
         if( _rationalLP != 0 )
            _rationalLP->changeSense(value == SoPlex2::OBJSENSE_MAXIMIZE ? SPxLPRational::MAXIMIZE : SPxLPRational::MINIMIZE);
         _invalidateSolution();
         break;

      // type of computational form, i.e., column or row representation
      case SoPlex2::REPRESENTATION:
         if( value != SoPlex2::REPRESENTATION_COLUMN && value != SoPlex2::REPRESENTATION_ROW && value != SoPlex2::REPRESENTATION_AUTO )
            return false;
         break;

      // type of algorithm, i.e., enter or leave
      case SoPlex2::ALGORITHM:
         if( value != SoPlex2::ALGORITHM_ENTER && value != SoPlex2::ALGORITHM_LEAVE )
            return false;
         _solver.setType(value == SoPlex2::ALGORITHM_ENTER ? SPxSolver::ENTER : SPxSolver::LEAVE);
         break;

      // type of LU update
      case SoPlex2::FACTOR_UPDATE_TYPE:
         if( value != SoPlex2::FACTOR_UPDATE_TYPE_ETA || value != SoPlex2::FACTOR_UPDATE_TYPE_FT )
            return false;
         _slufactor.setUtype(value == SoPlex2::FACTOR_UPDATE_TYPE_ETA ? SLUFactor::ETA : SLUFactor::FOREST_TOMLIN);
         break;

      // maximum number of updates before fresh factorization
      case SoPlex2::FACTOR_UPDATE_MAX:
         if( value <= 1 )
            return false;
         else
         {
            ///@todo set value in factorization
            break;
         }

      // iteration limit (-1 if unlimited)
      case SoPlex2::ITERLIMIT:
         if( value < -1 )
            return false;
         else
            break;

      // display frequency
      case SoPlex2::DISPLAYFREQ:
         if( value <= 0 )
            return false;
         else
         {
            ///@todo set value in solver
            break;
         }

      // verbosity level
      case SoPlex2::VERBOSITY:
         Param::setVerbose(value);
         break;

      // type of simplifier
      case SoPlex2::SIMPLIFIER:
         switch( value )
         {
         case SIMPLIFIER_OFF:
            _simplifier = 0;
            break;
         case SIMPLIFIER_AUTO:
            _simplifier = &_simplifierMainSM;
            assert(_simplifier != 0);
            break;
         default:
            return false;
         }
         break;

      // type of scaler
      case SoPlex2::SCALER:
         switch( value )
         {
         case SCALER_OFF:
            _scaler = 0;
            break;
         case SCALER_UNIEQUI:
            _scaler = &_scalerUniequi;
            break;
         case SCALER_BIEQUI:
            _scaler = &_scalerBiequi;
            break;
         case SCALER_GEO1:
            _scaler = &_scalerGeo1;
            break;
         case SCALER_GEO8:
            _scaler = &_scalerGeo8;
            break;
         default:
            return false;
         }
         break;

      // type of starter used to create crash basis
      case SoPlex2::STARTER:
         switch( value )
         {
         case STARTER_OFF:
            _starter = 0;
            break;
         case STARTER_WEIGHT:
            _starter = &_starterWeight;
            break;
         case STARTER_SUM:
            _starter = &_starterSum;
            break;
         case STARTER_VECTOR:
            _starter = &_starterVector;
            break;
         default:
            return false;
         }
         break;

      // type of pricer
      case SoPlex2::PRICER:
         switch( value )
         {
         case PRICER_AUTO:
            ///@todo implement in the solve routine
            _solver.setPricer(&_pricerQuickSteep);
            break;
         case PRICER_DANTZIG:
            _solver.setPricer(&_pricerDantzig);
            break;
         case PRICER_PARMULT:
            _solver.setPricer(&_pricerParMult);
            break;
         case PRICER_DEVEX:
            _solver.setPricer(&_pricerDevex);
            break;
         case PRICER_QUICKSTEEP:
            _solver.setPricer(&_pricerQuickSteep);
            break;
         case PRICER_STEEP:
            _solver.setPricer(&_pricerSteep);
            break;
         case PRICER_HYBRID:
            _solver.setPricer(&_pricerHybrid);
            break;
         default:
            return false;
         }
         break;

      // mode for synchronizing real and rational LP
      case SoPlex2::SYNCMODE:
         switch( value )
         {
         case SYNCMODE_ONLYREAL:
            if( _rationalLP != 0 )
            {
               _rationalLP->~SPxLPRational();
               spx_free(_rationalLP);
            }
            break;
         case SYNCMODE_AUTO:
            if( intParam(param) == SYNCMODE_ONLYREAL )
               _syncLPRational();
            break;
         case SYNCMODE_MANUAL:
            _ensureRationalLP();
            break;
         default:
            return false;
         }
         break;

      // mode for reading LP files; nothing to do but change the value if valid
      case SoPlex2::READMODE:
         switch( value )
         {
         case READMODE_REAL:
         case READMODE_RATIONAL:
            break;
         default:
            return false;
         }
         break;

      // mode for iterative refinement strategy; nothing to do but change the value if valid
      case SoPlex2::SOLVEMODE:
         switch( value )
         {
         case SOLVEMODE_REAL:
         case SOLVEMODE_AUTO:
         case SOLVEMODE_RATIONAL:
            break;
         default:
            return false;
         }
         break;

      // type of ratio test
      case SoPlex2::RATIOTESTER:
         switch( value )
         {
         case RATIOTESTER_TEXTBOOK:
            _solver.setTester(&_ratiotesterTextbook);
            break;
         case RATIOTESTER_HARRIS:
            _solver.setTester(&_ratiotesterHarris);
            break;
         case RATIOTESTER_FAST:
            _solver.setTester(&_ratiotesterFast);
            break;
         case RATIOTESTER_BOUNDFLIPPING:
            _solver.setTester(&_ratiotesterBoundFlipping);
            break;
         default:
            return false;
         }
         break;

      default:
         return false;
      }

      _currentSettings->_intParamValues[param] = value;
      return true;
   }



   /// sets real parameter value; returns true on success
   bool SoPlex2::setRealParam(const RealParam param, const Real value, const bool quiet, const bool init)
   {
      assert(param >= 0);
      assert(param < REALPARAM_COUNT);
      assert(init || _isConsistent());

      if( !init && value == realParam(param) )
         return true;

      if( value < _currentSettings->_realParamLower[param] || value > _currentSettings->_realParamUpper[param] )
         return false;

      switch( param )
      {
      // general zero tolerance
      case SoPlex2::EPSILON_ZERO:
         Param::setEpsilon(value);
         break;

      // zero tolerance used in factorization
      case SoPlex2::EPSILON_FACTORIZATION:
         Param::setEpsilonFactorization(value);
         break;

      // zero tolerance used in update of the factorization
      case SoPlex2::EPSILON_UPDATE:
         Param::setEpsilonUpdate(value);
         break;

      // infinity threshold
      case SoPlex2::INFTY:
         break;

      // time limit in seconds (INFTY if unlimited)
      case SoPlex2::TIMELIMIT:
         break;

      // lower limit on objective value is set in solveReal()
      case SoPlex2::OBJLIMIT_LOWER:
         break;

      // upper limit on objective value is set in solveReal()
      case SoPlex2::OBJLIMIT_UPPER:
         break;

      // working tolerance for feasibility in floating-point solver
      case SoPlex2::FPFEASTOL:
         break;

      // working tolerance for optimality in floating-point solver
      case SoPlex2::FPOPTTOL:
         break;

      // maximum increase of scaling factors between refinements
      case SoPlex2::MAXSCALEINCR:
         break;

      // lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)
      case SoPlex2::LIFTMINVAL:
         break;

      // upper threshold in lifting (nonzero matrix coefficients with larger absolute value will be reformulated)
      case SoPlex2::LIFTMAXVAL:
         break;

      default:
         return false;
      }

      _currentSettings->_realParamValues[param] = value;
      return true;
   }



   /// sets rational parameter value; returns true on success
   bool SoPlex2::setRationalParam(const RationalParam param, const Rational value, const bool quiet, const bool init)
   {
      assert(param >= 0);
      assert(param < RATIONALPARAM_COUNT);
      assert(init || _isConsistent());

      if( !init && value == rationalParam(param) )
         return true;

      if( value < _currentSettings->_rationalParamLower[param] || value > _currentSettings->_rationalParamUpper[param] )
         return false;

      switch( param )
      {
      // primal feasibility tolerance
      case SoPlex2::FEASTOL:
         _solver.setFeastol(Real(value));
         break;

      // dual feasibility tolerance
      case SoPlex2::OPTTOL:
         _solver.setOpttol(Real(value));
         break;

      default:
         return false;
      }

      _currentSettings->_rationalParamValues[param] = value;
      return true;
   }



   /// sets parameter settings; returns true on success
   bool SoPlex2::setSettings(const Settings& settings, const bool quiet, const bool init)
   {
      assert(init || _isConsistent());

      bool success = true;

      *_currentSettings = settings;

      for( int i = 0; i < SoPlex2::BOOLPARAM_COUNT; i++ )
         success &= setBoolParam((BoolParam)i, _currentSettings->_boolParamValues[i], quiet, init);

      for( int i = 0; i < SoPlex2::INTPARAM_COUNT; i++ )
         success &= setIntParam((IntParam)i, _currentSettings->_intParamValues[i], quiet, init);

      for( int i = 0; i < SoPlex2::REALPARAM_COUNT; i++ )
         success &= setRealParam((RealParam)i, _currentSettings->_realParamValues[i], quiet, init);

      for( int i = 0; i < SoPlex2::RATIONALPARAM_COUNT; i++ )
         success &= setRationalParam((RationalParam)i, _currentSettings->_rationalParamValues[i], quiet, init);

      assert(_isConsistent());

      return success;
   }



   /// writes settings file; returns true on success
   bool SoPlex2::saveSettingsFile(const char* filename, const bool onlyChanged) const
   {
      assert(filename != 0);

      std::ofstream file(filename);
      if( file == 0 )
         return false;

      file.setf(std::ios::left);
      file << "# SoPlex version " << SOPLEX_VERSION / 100 << "." << (SOPLEX_VERSION / 10) % 10 << "." << SOPLEX_VERSION % 10 << "." << SOPLEX_SUBVERSION << "\n";

      for( int i = 0; i < SoPlex2::BOOLPARAM_COUNT; i++ )
      {
         if( onlyChanged && _currentSettings->_boolParamValues[i] == _currentSettings->_boolParamDefault[i] )
            continue;

         file << "\n";
         file << "# " << _currentSettings->_boolParamDescription[i] << "\n";
         file << "# range {true, false}, default " << (_currentSettings->_boolParamDefault[i] ? "true\n" : "false\n");
         file << "bool:" << _currentSettings->_boolParamName[i] << " = " << (_currentSettings->_boolParamValues[i] ? "true\n" : "false\n");
      }

      for( int i = 0; i < SoPlex2::INTPARAM_COUNT; i++ )
      {
         if( onlyChanged && _currentSettings->_intParamValues[i] == _currentSettings->_intParamDefault[i] )
            continue;

         file << "\n";
         file << "# " << _currentSettings->_intParamDescription[i] << "\n";
         file << "# range [-2147483648,2147483647], default " << _currentSettings->_intParamDefault[i] << "\n";
         file << "int:" << _currentSettings->_intParamName[i] << " = " << _currentSettings->_intParamValues[i] << "\n";
      }

      for( int i = 0; i < SoPlex2::REALPARAM_COUNT; i++ )
      {
         if( onlyChanged && _currentSettings->_realParamValues[i] == _currentSettings->_realParamDefault[i] )
            continue;

         file << "\n";
         file << "# " << _currentSettings->_realParamDescription[i] << "\n";
         file << "# range [" << _currentSettings->_realParamLower[i] << "," << _currentSettings->_realParamLower[i]
            << "], default " << _currentSettings->_realParamDefault[i] << "\n";
         file << "real:" << _currentSettings->_realParamName[i] << " = " << _currentSettings->_realParamValues[i] << "\n";
      }

      for( int i = 0; i < SoPlex2::RATIONALPARAM_COUNT; i++ )
      {
         if( onlyChanged && _currentSettings->_rationalParamValues[i] == _currentSettings->_rationalParamDefault[i] )
            continue;

         file << "\n";
         file << "# " << _currentSettings->_rationalParamDescription[i] << "\n";
         file << "# range [" << _currentSettings->_rationalParamLower[i] << "," << _currentSettings->_rationalParamLower[i]
            << "], default " << _currentSettings->_rationalParamDefault[i] << "\n";
         file << "rational:" << _currentSettings->_rationalParamName[i] << " = " << _currentSettings->_rationalParamValues[i] << "\n";
      }

      return true;
   }



   /// reads settings file; returns true on success
   bool SoPlex2::loadSettingsFile(const char* filename)
   {
      assert(filename != 0);

      // start timing
      _statistics->readingTime.start();

      MSG_INFO1( spxout << "Loading settings file <" << filename << "> . . .\n" );

      // open file
      spxifstream file(filename);

      if( !file )
      {
         MSG_ERROR( spxout << "Error opening settings file.\n" );
         return false;
      }

      // read file
      char line[SET_MAX_LINE_LEN];
      int lineNumber = 0;
      bool readError = false;
      bool parseError = false;

      while( !readError && !parseError)
      {
         lineNumber++;
         readError = !file.getline(line, sizeof(line));
         if( !readError )
            parseError = !_parseSettingsLine(line, lineNumber);
      }

      if( readError && strlen(line) == SET_MAX_LINE_LEN - 1 )
      {
         MSG_ERROR( spxout << "Error reading settings file: line " << lineNumber << " in settings file exceeds " << SET_MAX_LINE_LEN - 2 << " characters.\n" );
      }

      // stop timing
      _statistics->readingTime.stop();

      return !readError && !parseError;
   }

   /// parses one setting string and returns true on success; note that string is modified
   bool SoPlex2::parseSettingsString(char* line)
   {
      assert(line != 0);

      // find the start of the parameter type
      while( *line == ' ' || *line == '\t' || *line == '\r' )
         line++;
      if( *line == '\0' || *line == '\n' || *line == '#' )
         return true;
      char* paramTypeString = line;

      // find the end of the parameter type
      while( *line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#' && *line != '\0' && *line != ':' )
         line++;
      if( *line == ':' )
      {
         *line = '\0';
         line++;
      }
      else
      {
         *line = '\0';
         line++;

         // search for the ':' char in the line
         while( *line == ' ' || *line == '\t' || *line == '\r' )
            line++;
         if( *line != ':' )
         {
            MSG_ERROR( spxout << "Error parsing setting string: no ':' separating parameter type and name.\n" );
            return false;
         }
         line++;
      }

      // find the start of the parameter name
      while( *line == ' ' || *line == '\t' || *line == '\r' )
         line++;
      if( *line == '\0' || *line == '\n' || *line == '#' )
      {
         MSG_ERROR( spxout << "Error parsing setting string: no parameter name.\n");
         return false;
      }
      char* paramName = line;

      // find the end of the parameter name
      while( *line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#' && *line != '\0' && *line != '=' )
         line++;
      if( *line == '=' )
      {
         *line = '\0';
         line++;
      }
      else
      {
         *line = '\0';
         line++;

         // search for the '=' char in the line
         while( *line == ' ' || *line == '\t' || *line == '\r' )
            line++;
         if( *line != '=' )
         {
            MSG_ERROR( spxout << "Error parsing setting string: no '=' after parameter name.\n" );
            return false;
         }
         line++;
      }

      // find the start of the parameter value string
      while( *line == ' ' || *line == '\t' || *line == '\r' )
         line++;
      if( *line == '\0' || *line == '\n' || *line == '#' )
      {
         MSG_ERROR( spxout << "Error parsing setting string: no parameter value.\n");
         return false;
      }
      char* paramValueString = line;

      // find the end of the parameter value string
      while( *line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#' && *line != '\0' )
         line++;
      if( *line != '\0' )
      {
         // check, if the rest of the line is clean
         *line = '\0';
         line++;
         while( *line == ' ' || *line == '\t' || *line == '\r' )
            line++;
         if( *line != '\0' && *line != '\n' && *line != '#' )
         {
            MSG_ERROR( spxout << "Error parsing setting string: additional character '" << *line << "' after parameter value.\n" );
            return false;
         }
      }

      // check whether we have a bool parameter
      if( strncmp(paramTypeString, "bool", 4) == 0 )
      {
         for( int param = 0; ; param++ )
         {
            if( param >= SoPlex2::BOOLPARAM_COUNT )
            {
               MSG_ERROR( spxout << "Error parsing setting string: unknown parameter name <" << paramName << ">.\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_boolParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               if( strncasecmp(paramValueString, "true", 4) == 0 )
               {
                  setBoolParam((SoPlex2::BoolParam)param, true);
                  break;
               }
               else if( strncasecmp(paramValueString, "false", 5) == 0 )
               {
                  setBoolParam((SoPlex2::BoolParam)param, false);
                  break;
               }
               else
               {
                  MSG_ERROR( spxout << "Error parsing setting string: invalid value <" << paramValueString << "> for bool parameter <" << paramName << ">.\n" );
                  return false;
               }
            }
         }

         return true;
      }

      // check whether we have an integer parameter
      if( strncmp(paramTypeString, "int", 3) == 0 )
      {
         for( int param = 0; ; param++ )
         {
            if( param >= SoPlex2::INTPARAM_COUNT )
            {
               MSG_ERROR( spxout << "Error parsing setting string: unknown parameter name <" << paramName << ">.\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_intParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               int value;

               if( sscanf(paramValueString, "%d", &value) == 1 && setIntParam((SoPlex2::IntParam)param, value) )
                  break;
               else
               {
                  MSG_ERROR( spxout << "Error parsing setting string: invalid value <" << paramValueString << "> for int parameter <" << paramName << ">.\n" );
                  return false;
               }
            }
         }

         return true;
      }

      // check whether we have a real parameter
      if( strncmp(paramTypeString, "real", 4) == 0 )
      {
         for( int param = 0; ; param++ )
         {
            if( param >= SoPlex2::REALPARAM_COUNT )
            {
               MSG_ERROR( spxout << "Error parsing setting string: unknown parameter name <" << paramName << ">.\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_realParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               Real value;

               if( sscanf(paramValueString, "%" REAL_FORMAT, &value) == 1 && setRealParam((SoPlex2::RealParam)param, value) )
                  break;
               else
               {
                  MSG_ERROR( spxout << "Error parsing setting string: invalid value <" << paramValueString << "> for real parameter <" << paramName << ">.\n" );
                  return false;
               }
            }
         }

         return true;
      }

      // check whether we have a rational parameter
      if( strncmp(paramTypeString, "rational", 8) == 0 )
      {
         for( int param = 0; ; param++ )
         {
            if( param >= SoPlex2::RATIONALPARAM_COUNT )
            {
               MSG_ERROR( spxout << "Error parsing setting string: unknown parameter name <" << paramName << ">.\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_rationalParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               Rational value;

               if( readStringRational(paramValueString, value) && setRationalParam((SoPlex2::RationalParam)param, value) )
                  break;
               else
               {
                  MSG_ERROR( spxout << "Error parsing setting string: invalid value <" << paramValueString << "> for rational parameter <" << paramName << ">.\n" );
                  return false;
               }
            }
         }

         return true;
      }

      MSG_ERROR( spxout << "Error parsing setting string: invalid parameter type <" << paramTypeString << "> for parameter <" << paramName << ">.\n" );

      return false;
   }




   /// prints statistics on real solution
   void SoPlex2::printSolutionStatisticsReal(std::ostream& os)
   {
      os << std::scientific << std::setprecision(8)
         << "Solution           : \n"
         << "  Value            : " << objValueReal() << "\n"
         << "  Proven primal    : " << "?\n"
         << "  Proven dual      : " << "?\n";
   }



   /// prints statistics on rational solution
   void SoPlex2::printSolutionStatisticsRational(std::ostream& os)
   {
      os << "Solution           : \n"
         << "  Objective value  : " << rationalToString(objValueRational()) << "\n"
         << "  Proven primal    : " << "?\n"
         << "  Proven dual      : " << "?\n";

      os << "Violations         : \n"
         << "  Max. bound       : " << "?\n"
         << "  Max. row         : " << "?\n"
         << "  Max. dual        : " << "?\n";
   }



   /// prints statistics on solving process
   void SoPlex2::printSolvingStatistics(std::ostream& os)
   {
      assert(_statistics != 0);
      _statistics->print(os);
   }



   /// prints complete real statistics
   void SoPlex2::printStatisticsReal(std::ostream& os)
   {
      os << std::setprecision(2);

      printStatus(os, _status);

      os << "Real LP            : \n"
         << "  Objective sense  : " << (intParam(SoPlex2::OBJSENSE) == SoPlex2::OBJSENSE_MINIMIZE ? "minimize\n" : "maximize\n");
      _realLP->printProblemStatistics(os);

      printSolutionStatisticsReal(os);

      printSolvingStatistics(os);
   }



   /// prints complete rational statistics
   void SoPlex2::printStatisticsRational(std::ostream& os)
   {
      os << std::setprecision(2);

      printStatus(os, _status);

      os << "Rational LP        : \n"
         << "  Objective sense  : " << (intParam(SoPlex2::OBJSENSE) == SoPlex2::OBJSENSE_MINIMIZE ? "minimize\n" : "maximize\n");
      _rationalLP->printProblemStatistics(os);

      printSolutionStatisticsRational(os);

      printSolvingStatistics(os);
   }



   /// prints status
   void SoPlex2::printStatus(std::ostream& os, SPxSolver::Status status)
   {
      os << "SoPlex status      : ";

      switch( status )
      {
      case SPxSolver::ERROR:
         os << "error [unspecified]";
         break;
      case SPxSolver::NO_RATIOTESTER:
         os << "error [no ratiotester loaded]";
         break;
      case SPxSolver::NO_PRICER:
         os << "error [no pricer loaded]";
         break;
      case SPxSolver::NO_SOLVER:
         os << "error [no linear solver loaded]";
         break;
      case SPxSolver::NOT_INIT:
         os << "error [not initialized]";
         break;
      case SPxSolver::ABORT_CYCLING:
         os << "solving aborted [cycling]";
         break;
      case SPxSolver::ABORT_TIME:
         os << "solving aborted [time limit reached]";
         break;
      case SPxSolver::ABORT_ITER:
         os << "solving aborted [iteration limit reached]";
         break;
      case SPxSolver::ABORT_VALUE:
         os << "solving aborted [objective limit reached]";
         break;
      case SPxSolver::NO_PROBLEM:
         os << "no problem loaded";
         break;
      case SPxSolver::REGULAR:
         os << "basis is regular";
         break;
      case SPxSolver::SINGULAR:
         os << "basis is singular";
         break;
      case SPxSolver::OPTIMAL:
         os << "problem is solved [optimal]";
         break;
      case SPxSolver::UNBOUNDED:
         os << "problem is solved [unbounded]";
         break;
      case SPxSolver::INFEASIBLE:
         os << "problem is solved [infeasible]";
         break;
      case SPxSolver::INForUNBD:
         os << "problem is solved [infeasible or unbounded]";
         break;
      default:
      case SPxSolver::UNKNOWN:
         os << "unknown";
         break;
      }

      os << "\n";
   }


   /// checks if real LP and rational LP are in sync; dimensions will always be compared,
   /// vector and matrix values only if the respective parameter is set to true.
   /// If quiet is set to true the function will only display which vectors are different.
   bool SoPlex2::areLPsInSync(const bool checkVecVals, const bool checkMatVals, const bool quiet) const
   {
      bool result = true;
      bool nRowsMatch = true;
      bool nColsMatch = true;
      bool rhsDimMatch = true;
      bool lhsDimMatch = true;
      bool maxObjDimMatch = true;
      bool upperDimMatch = true;
      bool lowerDimMatch = true;

      // compare number of Rows
      if( _realLP->nRows() != _rationalLP->nRows() )
      {
         MSG_ERROR( spxout << "The number of Rows in the Real LP does not match the one in the Rational LP."
               << " Real LP: " << _realLP->nRows() << "  Rational LP: " << _rationalLP->nRows() << std::endl);
         result = false;
         nRowsMatch = false;
      }

      // compare number of Columns
      if( _realLP->nCols() != _rationalLP->nCols() )
      {
         MSG_ERROR( spxout << "The number of Columns in the Real LP does not match the one in the Rational LP."
               << " Real LP: " << _realLP->nCols() << "  Rational LP: " << _rationalLP->nCols() << std::endl);
         result = false;
         nColsMatch = false;
      }

      // compare number of nonZeros
      if( _realLP->nNzos() != _rationalLP->nNzos() )
      {
         MSG_ERROR( spxout << "The number of nonZeros in the Real LP does not match the one in the Rational LP."
               << " Real LP: " << _realLP->nNzos() << "  Rational LP: " << _rationalLP->nNzos() << std::endl);
         result = false;
      }

      // compare the dimensions of the right hand side vectors
      if( _realLP->rhs().dim() != _rationalLP->rhs().dim() )
      {
         MSG_ERROR( spxout << "The dimension of the right hand side vector of the Real LP does not match the one of the Rational LP."
               << " Real LP: " << _realLP->rhs().dim() << "  Rational LP: " << _rationalLP->rhs().dim() << std::endl);
         result = false;
         rhsDimMatch = false;

      }

      // compare the dimensions of the left hand side vectors
      if( _realLP->lhs().dim() != _rationalLP->lhs().dim() )
      {
         MSG_ERROR( spxout << "The dimension of the left hand side vector of the Real LP does not match the one of the Rational LP."
               << " Real LP: " << _realLP->lhs().dim() << "  Rational LP: " << _rationalLP->lhs().dim() << std::endl);
         result = false;
         lhsDimMatch = false;
      }

      // compare the dimensions of the objective function vectors
      if( _realLP->maxObj().dim() != _rationalLP->maxObj().dim() )
      {
         MSG_ERROR( spxout << "The dimension of the objective function vector of the Real LP does not match the one of the Rational LP."
               << " Real LP: " << _realLP->maxObj().dim() << "  Rational LP: " << _rationalLP->maxObj().dim() << std::endl);
         result = false;
         maxObjDimMatch = false;
      }

      // compare the sense
      if( (int)_realLP->spxSense() != (int)_rationalLP->spxSense() )
         {
            MSG_ERROR( spxout << "The objective function sense of the Real LP does not match the one of the Rational LP."
                  << " Real LP: " << _realLP->spxSense() << "  Rational LP: " << _rationalLP->spxSense() << std::endl);
            result = false;
         }

      // compare the dimensions of upper bound vectors
      if( _realLP->upper().dim() != _rationalLP->upper().dim() )
      {
         MSG_ERROR( spxout << "The dimension of the upper bound vector of the Real LP does not match the one of the Rational LP."
               << " Real LP: " << _realLP->upper().dim() << "  Rational LP: " << _rationalLP->upper().dim() << std::endl);
         result = false;
         upperDimMatch = false;
      }

      // compare the dimensions of the objective function vectors
      if( _realLP->lower().dim() != _rationalLP->lower().dim() )
      {
         MSG_ERROR( spxout << "The dimension of the lower bound vector of the Real LP does not match the one of the Rational LP."
               << " Real LP: " << _realLP->lower().dim() << "  Rational LP: " << _rationalLP->lower().dim() << std::endl);
         result = false;
         lowerDimMatch = false;
      }

      // compares the values of the rhs, lhs, maxObj, upper, lower vectors
      if( checkVecVals )
      {
         bool rhsValMatch = true;
         bool lhsValMatch = true;
         bool maxObjValMatch = true;
         bool upperValMatch = true;
         bool lowerValMatch = true;

         // compares the values of the right hand side vectors
         if( rhsDimMatch )
         {
            for( int i = 0; i < _realLP->rhs().dim(); i++ )
            {
               if( !_rationalLP->rhs()[i].isAdjacentTo(_realLP->rhs()[i]) )
               {
                  if( !quiet )
                  {
                     MSG_ERROR( spxout << "Entries number " << i << " of the right hand side vectors don't match."
                           << " Real LP: " << _realLP->rhs()[i] << "  Rational LP: " << _rationalLP->rhs()[i] << std::endl);
                  }
                  rhsValMatch = false;
                  result = false;
               }
            }

            if( !rhsValMatch && quiet )
            {
               MSG_ERROR( spxout << "The values of the right hand side vectors don't match." << std::endl );
            }
         }

         // compares the values of the left hand side vectors
         if( lhsDimMatch )
         {
            for( int i = 0; i < _realLP->lhs().dim(); i++ )
            {
               if( !_rationalLP->lhs()[i].isAdjacentTo(_realLP->lhs()[i]) )
               {
                  if( !quiet )
                  {
                     MSG_ERROR( spxout << "Entries number " << i << " of the left hand side vectors don't match."
                           << " Real LP: " << _realLP->lhs()[i] << "  Rational LP: " << _rationalLP->lhs()[i] << std::endl);
                  }
                  lhsValMatch = false;
                  result = false;
               }
            }

            if( !lhsValMatch && quiet )
            {
               MSG_ERROR( spxout << "The values of the left hand side vectors don't match." << std::endl );
            }
         }

         // compares the values of the objective function vectors
         if( maxObjDimMatch )
         {
            for( int i = 0; i < _realLP->maxObj().dim(); i++ )
            {
               if( !_rationalLP->maxObj()[i].isAdjacentTo(_realLP->maxObj()[i]) )
               {
                  if( !quiet )
                  {
                     MSG_ERROR( spxout << "Entries number " << i << " of the objective function vectors don't match."
                           << " Real LP: " << _realLP->maxObj()[i] << "  Rational LP: " << _rationalLP->maxObj()[i] << std::endl);
                  }
                  maxObjValMatch = false;
                  result = false;
               }
            }

            if( !maxObjValMatch && quiet )
            {
               MSG_ERROR( spxout << "The values of the objective function vectors don't match." << std::endl );
            }
         }

         // compares the values of the upper bound vectors
         if( upperDimMatch )
         {
            for( int i = 0; i < _realLP->upper().dim(); i++ )
            {
               if( !_rationalLP->upper()[i].isAdjacentTo(_realLP->upper()[i]) )
               {
                  if( !quiet )
                  {
                     MSG_ERROR( spxout << "Entries number " << i << " of the upper bound vectors don't match."
                           << " Real LP: " << _realLP->upper()[i] << "  Rational LP: " << _rationalLP->upper()[i] << std::endl);
                  }
                  upperValMatch = false;
                  result = false;
               }
            }

            if( !upperValMatch && quiet )
            {
               MSG_ERROR( spxout << "The values of the upper bound vectors don't match." << std::endl );
            }
         }

         // compares the values of the lower bound vectors
         if( lowerDimMatch )
         {
            for( int i = 0; i < _realLP->lower().dim(); i++ )
            {
               if( !_rationalLP->lower()[i].isAdjacentTo(_realLP->lower()[i]) )
               {
                  if( !quiet )
                  {
                     MSG_ERROR( spxout << "Entries number " << i << " of the lower bound vectors don't match."
                           << " Real LP: " << _realLP->lower()[i] << "  Rational LP: " << _rationalLP->lower()[i] << std::endl);
                  }
                  lowerValMatch = false;
                  result = false;
               }
            }

            if( !lowerValMatch && quiet )
            {
               MSG_ERROR( spxout << "The values of the lower bound vectors don't match." << std::endl );
            }
         }
      }

      // compare the values of the matrix
      if( checkMatVals && nRowsMatch && nColsMatch )
      {
         bool matrixValMatch = true;

         for( int i = 0; i < _realLP->nCols() ; i++ )
         {
            for( int j = 0;j < _realLP->nRows() ; j++ )
            {
               if( !_rationalLP->colVector(i)[j].isAdjacentTo(_realLP->colVector(i)[j]) )
               {
                  if( !quiet )
                  {
                     MSG_ERROR( spxout << "Entries number " << j << " of column number " << i << " don't match."
                           << " Real LP: " << _realLP->colVector(i)[j] << "  Rational LP: " << _rationalLP->colVector(i)[j] << std::endl);
                  }
                  matrixValMatch = false;
                  result = false;
               }
            }
         }

         if( !matrixValMatch && quiet )
         {
            MSG_ERROR( spxout << "The values of the matrices don't match." << std::endl );
         }
      }

      return result;
   }



   /// creates a permutation for removing rows/columns from an array of indices
   void SoPlex2::_idxToPerm(int* idx, int idxSize, int* perm, int permSize) const
   {
      assert(idx != 0);
      assert(idxSize >= 0);
      assert(perm != 0);
      assert(permSize >= 0);

      for( int i = 0; i < permSize; i++ )
         perm[i] = i;

      for( int i = 0; i < idxSize; i++ )
      {
         assert(idx[i] >= 0);
         assert(idx[i] < permSize);
         perm[idx[i]] = -1;
      }
   }



   /// creates a permutation for removing rows/columns from a range of indices
   void SoPlex2::_rangeToPerm(int start, int end, int* perm, int permSize) const
   {
      assert(perm != 0);
      assert(permSize >= 0);

      for( int i = 0; i < permSize; i++ )
         perm[i] = (i < start || i > end) ? i : -1;
   }



   /// checks consistency
   bool SoPlex2::_isConsistent() const
   {
      assert(_statistics != 0);
      assert(_currentSettings != 0);

      assert(_realLP != 0);
      assert(_rationalLP != 0 || intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL);

      assert(_realLP != &_solver || _isRealLPLoaded);
      assert(_realLP == &_solver || !_isRealLPLoaded);

      assert(!_hasBasis || _isRealLPLoaded || _basisStatusRows.size() == numRowsReal());
      assert(!_hasBasis || _isRealLPLoaded || _basisStatusCols.size() == numColsReal());

      return true;
   }



   /// should solving process be stopped?
   bool SoPlex2::_isSolveStopped() const
   {
      assert(_statistics != 0);

      return _statistics->solvingTime.userTime() >= realParam(TIMELIMIT)
         || (intParam(ITERLIMIT) >= 0 && _statistics->iterations >= intParam(ITERLIMIT))
         || (intParam(REFLIMIT) >= 0 && _statistics->refinements >= intParam(REFLIMIT))
         || (intParam(STALLREFLIMIT) >= 0 && _statistics->stallRefinements >= intParam(STALLREFLIMIT));
   }



   /// parses one line in a settings file and returns true on success; note that the string is modified
   bool SoPlex2::_parseSettingsLine(char* line, const int lineNumber)
   {
      assert(line != 0);

      // find the start of the parameter type
      while( *line == ' ' || *line == '\t' || *line == '\r' )
         line++;
      if( *line == '\0' || *line == '\n' || *line == '#' )
         return true;
      char* paramTypeString = line;

      // find the end of the parameter type
      while( *line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#' && *line != '\0' && *line != ':' )
         line++;
      if( *line == ':' )
      {
         *line = '\0';
         line++;
      }
      else
      {
         *line = '\0';
         line++;

         // search for the ':' char in the line
         while( *line == ' ' || *line == '\t' || *line == '\r' )
            line++;
         if( *line != ':' )
         {
            MSG_ERROR( spxout << "Error parsing settings file: no ':' separating parameter type and name in line " << lineNumber << ".\n" );
            return false;
         }
         line++;
      }

      // find the start of the parameter name
      while( *line == ' ' || *line == '\t' || *line == '\r' )
         line++;
      if( *line == '\0' || *line == '\n' || *line == '#' )
      {
         MSG_ERROR( spxout << "Error parsing settings file: no parameter name in line " << lineNumber << ".\n");
         return false;
      }
      char* paramName = line;

      // find the end of the parameter name
      while( *line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#' && *line != '\0' && *line != '=' )
         line++;
      if( *line == '=' )
      {
         *line = '\0';
         line++;
      }
      else
      {
         *line = '\0';
         line++;

         // search for the '=' char in the line
         while( *line == ' ' || *line == '\t' || *line == '\r' )
            line++;
         if( *line != '=' )
         {
            MSG_ERROR( spxout << "Error parsing settings file: no '=' after parameter name in line " << lineNumber << ".\n" );
            return false;
         }
         line++;
      }

      // find the start of the parameter value string
      while( *line == ' ' || *line == '\t' || *line == '\r' )
         line++;
      if( *line == '\0' || *line == '\n' || *line == '#' )
      {
         MSG_ERROR( spxout << "Error parsing settings file: no parameter value in line " << lineNumber << ".\n");
         return false;
      }
      char* paramValueString = line;

      // find the end of the parameter value string
      while( *line != ' ' && *line != '\t' && *line != '\r' && *line != '\n' && *line != '#' && *line != '\0' )
         line++;
      if( *line != '\0' )
      {
         // check, if the rest of the line is clean
         *line = '\0';
         line++;
         while( *line == ' ' || *line == '\t' || *line == '\r' )
            line++;
         if( *line != '\0' && *line != '\n' && *line != '#' )
         {
            MSG_ERROR( spxout << "Error parsing settings file: additional character '" << *line << "' after parameter value in line " << lineNumber << ".\n" );
            return false;
         }
      }

      // check whether we have a bool parameter
      if( strncmp(paramTypeString, "bool", 4) == 0 )
      {
         for( int param = 0; ; param++ )
         {
            if( param >= SoPlex2::BOOLPARAM_COUNT )
            {
               MSG_ERROR( spxout << "Error parsing settings file: unknown parameter name <" << paramName << "> in line " << lineNumber << ".\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_boolParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               if( strncasecmp(paramValueString, "true", 4) == 0 )
               {
                  setBoolParam((SoPlex2::BoolParam)param, true);
                  break;
               }
               else if( strncasecmp(paramValueString, "false", 5) == 0 )
               {
                  setBoolParam((SoPlex2::BoolParam)param, false);
                  break;
               }
               else
               {
                  MSG_ERROR( spxout << "Error parsing settings file: invalid value <" << paramValueString << "> for bool parameter <" << paramName << "> in line " << lineNumber << ".\n" );
                  return false;
               }
            }
         }

         return true;
      }

      // check whether we have an integer parameter
      if( strncmp(paramTypeString, "int", 3) == 0 )
      {
         for( int param = 0; ; param++ )
         {
            if( param >= SoPlex2::INTPARAM_COUNT )
            {
               MSG_ERROR( spxout << "Error parsing settings file: unknown parameter name <" << paramName << "> in line " << lineNumber << ".\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_intParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               int value;

               if( sscanf(paramValueString, "%d", &value) == 1 && setIntParam((SoPlex2::IntParam)param, value) )
                  break;
               else
               {
                  MSG_ERROR( spxout << "Error parsing settings file: invalid value <" << paramValueString << "> for int parameter <" << paramName << "> in line " << lineNumber << ".\n" );
                  return false;
               }
            }
         }

         return true;
      }

      // check whether we have a real parameter
      if( strncmp(paramTypeString, "real", 4) == 0 )
      {
         for( int param = 0; ; param++ )
         {
            if( param >= SoPlex2::REALPARAM_COUNT )
            {
               MSG_ERROR( spxout << "Error parsing settings file: unknown parameter name <" << paramName << "> in line " << lineNumber << ".\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_realParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               Real value;

               if( sscanf(paramValueString, "%" REAL_FORMAT, &value) == 1 && setRealParam((SoPlex2::RealParam)param, value) )
                  break;
               else
               {
                  MSG_ERROR( spxout << "Error parsing settings file: invalid value <" << paramValueString << "> for real parameter <" << paramName << "> in line " << lineNumber << ".\n" );
                  return false;
               }
            }
         }

         return true;
      }

      // check whether we have a rational parameter
      if( strncmp(paramTypeString, "rational", 8) == 0 )
      {
         for( int param = 0; ; param++ )
         {
            if( param >= SoPlex2::RATIONALPARAM_COUNT )
            {
               MSG_ERROR( spxout << "Error parsing settings file: unknown parameter name <" << paramName << "> in line " << lineNumber << ".\n" );
               return false;
            }
            else if( strncmp(paramName, _currentSettings->_rationalParamName[param].c_str(), SET_MAX_LINE_LEN) == 0 )
            {
               Rational value;

               if( readStringRational(paramValueString, value) && setRationalParam((SoPlex2::RationalParam)param, value) )
                  break;
               else
               {
                  MSG_ERROR( spxout << "Error parsing settings file: invalid value <" << paramValueString << "> for rational parameter <" << paramName << "> in line " << lineNumber << ".\n" );
                  return false;
               }
            }
         }

         return true;
      }

      MSG_ERROR( spxout << "Error parsing settings file: invalid parameter type <" << paramTypeString << "> for parameter <" << paramName << "> in line " << lineNumber << ".\n" );

      return false;
   }



   /// adds a single row to the real LP and adjusts basis
   void SoPlex2::_addRowReal(const LPRowReal& lprow)
   {
      assert(_realLP != 0);

      _realLP->addRow(lprow);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
         _basisStatusRows.append(SPxSolver::BASIC);
   }



   /// adds multiple rows to the real LP and adjusts basis
   void SoPlex2::_addRowsReal(const LPRowSetReal& lprowset)
   {
      assert(_realLP != 0);

      _realLP->addRows(lprowset);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
         _basisStatusRows.append(lprowset.num(), SPxSolver::BASIC);
   }


   /// adds a single column to the real LP and adjusts basis
   void SoPlex2::_addColReal(const LPColReal& lpcol)
   {
      assert(_realLP != 0);

      _realLP->addCol(lpcol);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
      {
         if( lpcol.lower() > -realParam(SoPlex2::INFTY) )
            _basisStatusCols.append(SPxSolver::ON_LOWER);
         else if( lpcol.upper() < realParam(SoPlex2::INFTY) )
            _basisStatusCols.append(SPxSolver::ON_UPPER);
         else
            _basisStatusCols.append(SPxSolver::ZERO);
      }
   }



   /// adds multiple columns to the real LP and adjusts basis
   void SoPlex2::_addColsReal(const LPColSetReal& lpcolset)
   {
      assert(_realLP != 0);

      _realLP->addCols(lpcolset);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
      {
         for( int i = 0; i < lpcolset.num(); i++ )
         {
            if( lpcolset.lower(i) > -realParam(SoPlex2::INFTY) )
               _basisStatusCols.append(SPxSolver::ON_LOWER);
            else if( lpcolset.upper(i) < realParam(SoPlex2::INFTY) )
               _basisStatusCols.append(SPxSolver::ON_UPPER);
            else
               _basisStatusCols.append(SPxSolver::ZERO);
         }
      }
   }


   /// replaces row \p i with \p lprow and adjusts basis
   void SoPlex2::_changeRowReal(int i, const LPRowReal& lprow)
   {
      assert(_realLP != 0);

      _realLP->changeRow(i, lprow);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
      {
         if( _basisStatusRows[i] != SPxSolver::BASIC )
            _hasBasis = false;
         else if( _basisStatusRows[i] == SPxSolver::ON_LOWER && lprow.lhs() <= -realParam(SoPlex2::INFTY) )
            _basisStatusRows[i] = (lprow.rhs() < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusRows[i] == SPxSolver::ON_UPPER && lprow.rhs() >= realParam(SoPlex2::INFTY) )
            _basisStatusRows[i] = (lprow.lhs() > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }
   }



   /// changes left-hand side vector for constraints to \p lhs and adjusts basis
   void SoPlex2::_changeLhsReal(const VectorReal& lhs)
   {
      assert(_realLP != 0);

      _realLP->changeLhs(lhs);

      if( _isRealLPLoaded )
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      else if( _hasBasis )
      {
         for( int i = numRowsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusRows[i] == SPxSolver::ON_LOWER && lhs[i] <= -realParam(SoPlex2::INFTY) )
               _basisStatusRows[i] = (rhsReal(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         }
      }
   }



   /// changes left-hand side of row \p i to \p lhs and adjusts basis
   void SoPlex2::_changeLhsReal(int i, Real lhs)
   {
      assert(_realLP != 0);

      _realLP->changeLhs(i, lhs);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis && _basisStatusRows[i] == SPxSolver::ON_LOWER && lhs <= -realParam(SoPlex2::INFTY) )
         _basisStatusRows[i] = (rhsReal(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;

   }



   /// changes right-hand side vector to \p rhs and adjusts basis
   void SoPlex2::_changeRhsReal(const VectorReal& rhs)
   {
      assert(_realLP != 0);

      _realLP->changeRhs(rhs);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numRowsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusRows[i] == SPxSolver::ON_UPPER && rhs[i] >= realParam(SoPlex2::INFTY) )
               _basisStatusRows[i] = (lhsReal(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }
   }



   /// changes right-hand side of row \p i to \p rhs and adjusts basis
   void SoPlex2::_changeRhsReal(int i, Real rhs)
   {
      assert(_realLP != 0);

      _realLP->changeRhs(i, rhs);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis && _basisStatusRows[i] == SPxSolver::ON_UPPER && rhs >= realParam(SoPlex2::INFTY) )
         _basisStatusRows[i] = (lhsReal(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
   }



   /// changes left- and right-hand side vectors and adjusts basis
   void SoPlex2::_changeRangeReal(const VectorReal& lhs, const VectorReal& rhs)
   {
      assert(_realLP != 0);

      _realLP->changeRange(lhs, rhs);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numRowsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusRows[i] == SPxSolver::ON_LOWER && lhs[i] <= -realParam(SoPlex2::INFTY) )
               _basisStatusRows[i] = (rhs[i] < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
            else if( _basisStatusRows[i] == SPxSolver::ON_UPPER && rhs[i] >= realParam(SoPlex2::INFTY) )
               _basisStatusRows[i] = (lhs[i] > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }
   }



   /// changes left- and right-hand side of row \p i and adjusts basis
   void SoPlex2::_changeRangeReal(int i, Real lhs, Real rhs)
   {
      assert(_realLP != 0);

      _realLP->changeRange(i, lhs, rhs);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         if( _basisStatusRows[i] == SPxSolver::ON_LOWER && lhs <= -realParam(SoPlex2::INFTY) )
            _basisStatusRows[i] = (rhs < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusRows[i] == SPxSolver::ON_UPPER && rhs >= realParam(SoPlex2::INFTY) )
            _basisStatusRows[i] = (lhs > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }
   }



   /// replaces column \p i with \p lpcol and adjusts basis
   void SoPlex2::_changeColReal(int i, const LPColReal& lpcol)
   {
      assert(_realLP != 0);

      _realLP->changeCol(i, lpcol);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         if( _basisStatusCols[i] == SPxSolver::BASIC )
            _hasBasis = false;
         else if( _basisStatusCols[i] == SPxSolver::ON_LOWER && lpcol.lower() <= -realParam(SoPlex2::INFTY) )
            _basisStatusCols[i] = (lpcol.upper() < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusCols[i] == SPxSolver::ON_UPPER && lpcol.upper() >= realParam(SoPlex2::INFTY) )
            _basisStatusCols[i] = (lpcol.lower() > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }
   }



   /// changes vector of lower bounds to \p lower and adjusts basis
   void SoPlex2::_changeLowerReal(const VectorReal& lower)
   {
      assert(_realLP != 0);

      _realLP->changeLower(lower);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numColsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusCols[i] == SPxSolver::ON_LOWER && lower[i] <= -realParam(SoPlex2::INFTY) )
               _basisStatusCols[i] = (upperReal(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         }
      }
   }



   /// changes lower bound of column i to \p lower and adjusts basis
   void SoPlex2::_changeLowerReal(int i, Real lower)
   {
      assert(_realLP != 0);

      _realLP->changeLower(i, lower);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis && _basisStatusCols[i] == SPxSolver::ON_LOWER && lower <= -realParam(SoPlex2::INFTY) )
         _basisStatusCols[i] = (upperReal(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
   }



   /// changes vector of upper bounds to \p upper and adjusts basis
   void SoPlex2::_changeUpperReal(const VectorReal& upper)
   {
      assert(_realLP != 0);

      _realLP->changeUpper(upper);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numColsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusCols[i] == SPxSolver::ON_UPPER && upper[i] >= realParam(SoPlex2::INFTY) )
               _basisStatusCols[i] = (lowerReal(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }
   }



   /// changes \p i 'th upper bound to \p upper and adjusts basis
   void SoPlex2::_changeUpperReal(int i, Real upper)
   {
      assert(_realLP != 0);

      _realLP->changeUpper(i, upper);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis &&  _basisStatusCols[i] == SPxSolver::ON_UPPER && upper >= realParam(SoPlex2::INFTY) )
         _basisStatusCols[i] = (lowerReal(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
   }



   /// changes vectors of column bounds to \p lower and \p upper and adjusts basis
   void SoPlex2::_changeBoundsReal(const VectorReal& lower, const VectorReal& upper)
   {
      assert(_realLP != 0);

      _realLP->changeBounds(lower, upper);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numColsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusCols[i] == SPxSolver::ON_LOWER && lower[i] <= -realParam(SoPlex2::INFTY) )
               _basisStatusCols[i] = (upper[i] < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
            else if( _basisStatusCols[i] == SPxSolver::ON_UPPER && upper[i] >= realParam(SoPlex2::INFTY) )
               _basisStatusCols[i] = (lower[i] > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }
   }



   /// changes bounds of column \p i to \p lower and \p upper and adjusts basis
   void SoPlex2::_changeBoundsReal(int i, Real lower, Real upper)
   {
      assert(_realLP != 0);

      _realLP->changeBounds(i, lower, upper);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         if( _basisStatusCols[i] == SPxSolver::ON_LOWER && lower <= -realParam(SoPlex2::INFTY) )
            _basisStatusCols[i] = (upper < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusCols[i] == SPxSolver::ON_UPPER && upper >= realParam(SoPlex2::INFTY) )
            _basisStatusCols[i] = (lower > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }
   }



   /// changes matrix entry in row \p i and column \p j to \p val and adjusts basis
   void SoPlex2::_changeElementReal(int i, int j, Real val)
   {
      assert(_realLP != 0);

      _realLP->changeElement(i, j, val);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         if( _basisStatusRows[i] != SPxSolver::BASIC && _basisStatusCols[i] == SPxSolver::BASIC )
            _hasBasis = false;
      }
   }



   /// removes row \p i and adjusts basis
   void SoPlex2::_removeRowReal(int i)
   {
      assert(_realLP != 0);

      _realLP->removeRow(i);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         if( _basisStatusRows[i] != SPxSolver::BASIC )
            _hasBasis = false;
         else
         {
            _basisStatusRows[i] = _basisStatusRows[_basisStatusRows.size() - 1];
            _basisStatusRows.removeLast();
         }
      }
   }



   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsReal()
   void SoPlex2::_removeRowsReal(int perm[])
   {
      assert(_realLP != 0);

      _realLP->removeRows(perm);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numRowsReal() - 1; i >= 0 && _hasBasis; i-- )
         {
            if( perm[i] < 0 && _basisStatusRows[i] != SPxSolver::BASIC )
               _hasBasis = false;
            else if( perm[i] >= 0 && perm[i] != i )
            {
               assert(perm[i] < numRowsReal());
               assert(perm[perm[i]] < 0);

               _basisStatusRows[perm[i]] = _basisStatusRows[i];
            }
         }

         if( _hasBasis )
            _basisStatusRows.reSize(numRowsReal());
      }
   }



   /// removes column i
   void SoPlex2::_removeColReal(int i)
   {
      assert(_realLP != 0);

      _realLP->removeCol(i);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         if( _basisStatusCols[i] == SPxSolver::BASIC )
            _hasBasis = false;
         else
         {
            _basisStatusCols[i] = _basisStatusCols[_basisStatusCols.size() - 1];
            _basisStatusCols.removeLast();
         }
      }
   }



   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsReal()
   void SoPlex2::_removeColsReal(int perm[])
   {
      assert(_realLP != 0);

      _realLP->removeCols(perm);

      if( _isRealLPLoaded )
      {
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }
      else if( _hasBasis )
      {
         for( int i = numColsReal() - 1; i >= 0 && _hasBasis; i-- )
         {
            if( perm[i] < 0 && _basisStatusCols[i] == SPxSolver::BASIC )
               _hasBasis = false;
            else if( perm[i] >= 0 && perm[i] != i )
            {
               assert(perm[i] < numColsReal());
               assert(perm[perm[i]] < 0);

               _basisStatusCols[perm[i]] = _basisStatusCols[i];
            }
         }

         if( _hasBasis )
            _basisStatusCols.reSize(numColsReal());
      }
   }



   /// invalidates solution
   void SoPlex2::_invalidateSolution()
   {
      ///@todo maybe this should be done individually at the places when this method is called
      _status = SPxSolver::UNKNOWN;

      _solReal.invalidate();
      _hasSolReal = false;

      _solRational.invalidate();
      _hasSolRational = false;
   }



   /// enables simplifier and scaler
   void SoPlex2::_enableSimplifierAndScaler()
   {
      // type of simplifier
      switch( intParam(SoPlex2::SIMPLIFIER) )
      {
      case SIMPLIFIER_OFF:
         _simplifier = 0;
         break;
      case SIMPLIFIER_AUTO:
         _simplifier = &_simplifierMainSM;
         assert(_simplifier != 0);
         break;
      default:
         break;
      }

      // type of scaler
      switch( intParam(SoPlex2::SCALER) )
      {
      case SCALER_OFF:
         _scaler = 0;
         break;
      case SCALER_UNIEQUI:
         _scaler = &_scalerUniequi;
         break;
      case SCALER_BIEQUI:
         _scaler = &_scalerBiequi;
         break;
      case SCALER_GEO1:
         _scaler = &_scalerGeo1;
         break;
      case SCALER_GEO8:
         _scaler = &_scalerGeo8;
         break;
      default:
         break;
      }
   }



   /// disables simplifier and scaler
   void SoPlex2::_disableSimplifierAndScaler()
   {
      _simplifier = 0;
      _scaler = 0;
   }



   /// ensures that the rational LP is available; performs no sync
   void SoPlex2::_ensureRationalLP()
   {
      if( _rationalLP == 0 )
      {
         spx_alloc(_rationalLP);
         _rationalLP = new (_rationalLP) SPxLPRational();
      }
   }



   /// ensures that the real LP and the basis are loaded in the solver; performs no sync
   void SoPlex2::_ensureRealLPLoaded()
   {
      if( !_isRealLPLoaded )
      {
         assert(_realLP != &_solver);

         _solver.loadLP(*_realLP);
         _realLP->~SPxLPReal();
         spx_free(_realLP);
         _realLP = &_solver;
         _isRealLPLoaded = true;

         if( _hasBasis )
         {
            ///@todo this should not fail even if the basis is invalid (wrong dimension or wrong number of basic
            ///      entries); fix either in SPxSolver or in SPxBasis
            assert(_basisStatusRows.size() == numRowsReal());
            assert(_basisStatusCols.size() == numColsReal());
            _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
            _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
         }
      }
   }



   /// call floating-point solver and update statistics on iterations etc.
   void SoPlex2::_solveRealLPAndRecordStatistics()
   {
      bool _hadBasis = _hasBasis;

      // set time and iteration limit
      if( intParam(SoPlex2::ITERLIMIT) >= 0 )
         _solver.setTerminationIter(intParam(SoPlex2::ITERLIMIT) - _statistics->iterations);
      if( realParam(SoPlex2::TIMELIMIT) < realParam(SoPlex2::INFTY) )
         _solver.setTerminationTime(realParam(SoPlex2::TIMELIMIT) - _statistics->solvingTime.userTime());

      // ensure that tolerances are not too small
      if( _solver.feastol() < 1e-12 )
         _solver.setFeastol(1e-12);
      if( _solver.opttol() < 1e-12 )
         _solver.setOpttol(1e-12);

      // set correct representation
      if( (intParam(SoPlex2::REPRESENTATION) == SoPlex2::REPRESENTATION_COLUMN
            || (intParam(SoPlex2::REPRESENTATION) == SoPlex2::REPRESENTATION_AUTO && _solver.nCols() > _solver.nRows()))
            && _solver.rep() != SPxSolver::COLUMN )
      {
         _solver.setRep(SPxSolver::COLUMN);
      }
      else if( (intParam(SoPlex2::REPRESENTATION) == SoPlex2::REPRESENTATION_ROW
            || (intParam(SoPlex2::REPRESENTATION) == SoPlex2::REPRESENTATION_AUTO && _solver.nCols() < _solver.nRows()))
            && _solver.rep() != SPxSolver::ROW )
      {
         _solver.setRep(SPxSolver::ROW);
      }

      // call floating-point solver and catch exceptions
      _statistics->simplexTime.start();
      try
      {
         _solver.solve();
      }
      catch( SPxException E )
      {
         MSG_ERROR( spxout << "Caught exception <" << E.what() << "> while solving real LP.\n" );
      }
      catch( ... )
      {
         MSG_ERROR( spxout << "Caught unknown exception while solving real LP.\n" );
         _status = SPxSolver::ERROR;
      }
      _statistics->simplexTime.stop();

      // record statistics
      _statistics->iterations += _solver.iterations();
      _statistics->iterationsPrimal += _solver.primalIterations();
      _statistics->iterationsFromBasis += _hadBasis ? _solver.iterations() : 0;
      _statistics->luFactorizationTime += _slufactor.getFactorTime();
      _slufactor.resetFactorTime();
      _statistics->luSolveTime += _slufactor.getSolveTime();
      _slufactor.resetSolveTime();
      _statistics->luFactorizations += _slufactor.getFactorCount();
      _statistics->luSolves += _slufactor.getSolveCount();
   }



   /// reads real LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool SoPlex2::_readFileReal(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
   {
      assert(_realLP != 0);

      // clear statistics
      _statistics->clearAllData();

      // update status
      _status = SPxSolver::UNKNOWN;
      _invalidateSolution();
      _hasBasis = false;

      // start timing
      _statistics->readingTime.start();

      // read
      bool success = _realLP->readFile(filename, rowNames, colNames, intVars);
      setIntParam(SoPlex2::OBJSENSE, (_realLP->spxSense() == SPxLPReal::MAXIMIZE ? SoPlex2::OBJSENSE_MAXIMIZE : SoPlex2::OBJSENSE_MINIMIZE), true, true);

      // stop timing
      _statistics->readingTime.stop();

      if( success )
      {
         // if sync mode is auto, we have to copy the (rounded) real LP to the rational LP; this is counted to sync time
         // and not to reading time
         if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
            _syncLPRational();
      }
      else
         clearLPReal();

      return success;
   }



   /// reads rational LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool SoPlex2::_readFileRational(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
   {
      assert(_rationalLP != 0);

      // clear statistics
      _statistics->clearAllData();

      // start timing
      _statistics->readingTime.start();

      // update status
      _status = SPxSolver::UNKNOWN;
      _invalidateSolution();
      _hasBasis = false;

      // read
      _ensureRationalLP();
      bool success = _rationalLP->readFile(filename, rowNames, colNames, intVars);
      setIntParam(SoPlex2::OBJSENSE, (_rationalLP->spxSense() == SPxLPRational::MAXIMIZE ? SoPlex2::OBJSENSE_MAXIMIZE : SoPlex2::OBJSENSE_MINIMIZE), true, true);

      // stop timing
      _statistics->readingTime.stop();

      if( success )
      {
         // if sync mode is auto, we have to copy the (rounded) real LP to the rational LP; this is counted to sync time
         // and not to reading time
         if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_AUTO )
            _syncLPReal();
         // if a rational LP file is read, but only the (rounded) real LP should be kept, we have to free the rational LP
         else if( intParam(SoPlex2::SYNCMODE) == SYNCMODE_ONLYREAL )
         {
            _syncLPReal();
            _rationalLP->~SPxLPRational();
            spx_free(_rationalLP);
         }
      }
      else
         clearLPRational();

      return success;
   }



   /// synchronizes real LP with rational LP, i.e., copies (rounded) rational LP into real LP, without looking at the sync mode
   void SoPlex2::_syncLPReal()
   {
      // start timing
      _statistics->syncTime.start();

      // copy LP
      if( _isRealLPLoaded )
         _solver.loadLP((SPxLPReal)(*_rationalLP));
      else
         *_realLP = *_rationalLP;

      ///@todo try loading old basis
      _hasBasis = false;

      // invalidate real solution
      _solReal.invalidate();
      _hasSolReal = false;

      // stop timing
      _statistics->syncTime.stop();
   }



   /// synchronizes rational LP with real LP, i.e., copies real LP to rational LP, without looking at the sync mode
   void SoPlex2::_syncLPRational()
   {
      // start timing
      _statistics->syncTime.start();

      // copy LP
      _ensureRationalLP();
      *_rationalLP = *_realLP;

      // invalidate rational solution
      _solRational.invalidate();
      _hasSolRational = false;

      // stop timing
      _statistics->syncTime.stop();
   }



   /// synchronizes real solution with rational solution, i.e., copies real solution to rational solution
   void SoPlex2::_syncRealSolution()
   {
      if( _hasSolRational && !_hasSolReal )
      {
         _solReal = _solRational;
         _hasSolReal = true;
      }
   }



   /// synchronizes rational solution with real solution, i.e., copies (rounded) rational solution to real solution
   void SoPlex2::_syncRationalSolution()
   {
      if( _hasSolReal && !_hasSolRational )
      {
         _solRational = _solReal;
         _hasSolRational = true;
      }
   }
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
