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

/// maximum length of lines in settings file
#define SET_MAX_LINE_LEN 500

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

            // objective sense
            _intParamName[SoPlex2::OBJSENSE] = "objsense";
            _intParamDescription[SoPlex2::OBJSENSE] = "objective sense (-1 - minimize, +1 - maximize)";
            _intParamDefault[SoPlex2::OBJSENSE] = SoPlex2::OBJSENSE_MAXIMIZE;

            // type of computational form, i.e., column or row representation
            _intParamName[SoPlex2::REPRESENTATION] = "representation";
            _intParamDescription[SoPlex2::REPRESENTATION] = "type of computational form (0 - column representation, 1 - row representation)";
            _intParamDefault[SoPlex2::REPRESENTATION] = SoPlex2::REPRESENTATION_COLUMN;

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

            // type of scaler applied before simplification
            _intParamName[SoPlex2::SCALER_BEFORE_SIMPLIFIER] = "scaler_before_simplifier";
            _intParamDescription[SoPlex2::SCALER_BEFORE_SIMPLIFIER] = "scaling before simplification (0 - off, 1 - uni-equilibrium, 2 - bi-equilibrium, 3 - geometric, 4 - iterated geometric)";
            _intParamDefault[SoPlex2::SCALER_BEFORE_SIMPLIFIER] = SoPlex2::SCALER_BIEQUI;

            // type of scaler applied after simplification
            _intParamName[SoPlex2::SCALER_AFTER_SIMPLIFIER] = "scaler_after_simplifier";
            _intParamDescription[SoPlex2::SCALER_AFTER_SIMPLIFIER] = "scaling after simplification (0 - off, 1 - uni-equilibrium, 2 - bi-equilibrium, 3 - geometric, 4 - iterated geometric)";
            _intParamDefault[SoPlex2::SCALER_AFTER_SIMPLIFIER] = SoPlex2::SCALER_OFF;

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
            _realParamDefault[SoPlex2::MAXSCALEINCR] = DEFAULT_INFINITY;

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
      : _statistics(0)
      , _currentSettings(0)
      , _scalerUniequi(false)
      , _scalerBiequi(true)
      , _scalerGeo1(1)
      , _scalerGeo8(8)
      , _simplifier(0)
      , _firstScaler(0)
      , _secondScaler(0)
      , _starter(0)
      , _statusReal(SPxSolver::UNKNOWN)
      , _hasBasisReal(false)
      , _hasPrimalReal(false)
      , _hasPrimalrayReal(false)
      , _hasDualReal(false)
      , _hasDualfarkasReal(false)
      , _statusRational(SPxSolver::UNKNOWN)
      , _hasBasisRational(false)
   {
      // give lu factorization to solver
      _solver.setSolver(&_slufactor);

      // the real LP is initially stored in the solver
      _realLP = &_solver;
      _isRealLPLoaded = true;

      // construct the rational LP
      spx_alloc(_rationalLP);
      _rationalLP = new (_rationalLP) SPxLPRational();

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
      assert(rhs._rationalLP != 0);

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

         // copy basis
         _basisStatusRowsReal = rhs._basisStatusRowsReal;
         _basisStatusColsReal = rhs._basisStatusColsReal;
         _basisStatusRowsRational = rhs._basisStatusRowsRational;
         _basisStatusColsRational = rhs._basisStatusColsRational;

         // initialize pointers for simplifier, scalers, and starter
         setIntParam(SoPlex2::SIMPLIFIER, intParam(SoPlex2::SIMPLIFIER), true, true);
         setIntParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER, intParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER), true, true);
         setIntParam(SoPlex2::SCALER_AFTER_SIMPLIFIER, intParam(SoPlex2::SCALER_AFTER_SIMPLIFIER), true, true);
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
         _rationalLP = 0;
         spx_alloc(_rationalLP);
         _rationalLP = new (_rationalLP) SPxLPRational(*rhs._rationalLP);

         // copy boolean flags
         _isRealLPLoaded = rhs._isRealLPLoaded;
         _hasBasisReal = rhs._hasBasisReal;
         _hasBasisRational = rhs._hasBasisRational;
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
      assert(_realLP != 0);
      assert(_rationalLP != 0);

      // free real LP if different from the LP in the solver
      if( _realLP != &_solver )
      {
         spx_free(_realLP);
      }

      // free rational LP
      spx_free(_rationalLP);
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



   /// returns row identifier for row \p i
   SPxRowId SoPlex2::rowIdReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->rId(i);
   }



   /// returns column identifier for column \p i
   SPxColId SoPlex2::colIdReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->cId(i);
   }



   /// returns index of the row with identifier \p id
   int SoPlex2::idxReal(const SPxRowId& id) const
   {
      assert(_realLP != 0);
      return _realLP->number(id);
   }



   /// returns index of the column with identifier \p id
   int SoPlex2::idxReal(const SPxColId& id) const
   {
      assert(_realLP != 0);
      return _realLP->number(id);
   }



   /// returns index of the row or column with identifier \p id
   int SoPlex2::idxReal(const SPxId& id) const
   {
      assert(_realLP != 0);
      return _realLP->number(id);
   }



   /// gets row \p i
   void SoPlex2::getRowReal(int i, LPRowReal& lprow) const
   {
      assert(_realLP != 0);
      _realLP->getRow(i, lprow);
   }



   /// gets row with identifier \p id
   void SoPlex2::getRowReal(const SPxRowId& id, LPRowReal& lprow) const
   {
      assert(_realLP != 0);
      _realLP->getRow(id, lprow);
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



   /// returns vector of row with identifier \p id
   const SVectorReal& SoPlex2::rowVectorReal(const SPxRowId& id) const
   {
      assert(_realLP != 0);
      return _realLP->rowVector(id);
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



   /// returns right-hand side of row with identifier \p id
   Real SoPlex2::rhsReal(const SPxRowId& id) const
   {
      assert(_realLP != 0);
      return _realLP->rhs(id);
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



   /// returns left-hand side of row with identifier \p id
   Real SoPlex2::lhsReal(const SPxRowId& id) const
   {
      assert(_realLP != 0);
      return _realLP->lhs(id);
   }



   /// returns inequality type of row \p i
   LPRowReal::Type SoPlex2::rowTypeReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->rowType(i);
   }



   /// returns inequality type of row with identifier \p id
   LPRowReal::Type SoPlex2::rowTypeReal(const SPxRowId& id) const
   {
      assert(_realLP != 0);
      return _realLP->rowType(id);
   }



   /// gets column \p i
   void SoPlex2::getColReal(int i, LPColReal& lpcol) const
   {
      assert(_realLP != 0);
      return _realLP->getCol(i, lpcol);
   }



   /// gets column with identifier \p id.
   void SoPlex2::getColReal(const SPxColId& id, LPColReal& lpcol) const
   {
      assert(_realLP != 0);
      return _realLP->getCol(id, lpcol);
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



   /// returns vector of column with identifier \p id
   const SVectorReal& SoPlex2::colVectorReal(const SPxColId& id) const
   {
      assert(_realLP != 0);
      return _realLP->colVector(id);
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



   /// returns upper bound of column with identifier \p id
   Real SoPlex2::upperReal(const SPxColId& id) const
   {
      assert(_realLP != 0);
      return _realLP->upper(id);
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



   /// returns lower bound of column with identifier \p id
   Real SoPlex2::lowerReal(const SPxColId& id) const
   {
      assert(_realLP != 0);
      return _realLP->lower(id);
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



   /// returns objective value of column with identifier \p id
   Real SoPlex2::objReal(const SPxColId& id) const
   {
      assert(_realLP != 0);
      return _realLP->obj(id);
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



   /// returns objective value of column with identifier \p id after transformation to a maximization problem; since
   /// this is how it is stored internally, this is generally faster
   Real SoPlex2::maxObjReal(const SPxColId& id) const
   {
      assert(_realLP != 0);
      return _realLP->maxObj(id);
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



   /// returns row identifier for row \p i
   SPxRowId SoPlex2::rowIdRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->rId(i);
   }



   /// returns column identifier for column \p i
   SPxColId SoPlex2::colIdRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->cId(i);
   }



   /// returns index of the row with identifier \p id
   int SoPlex2::idxRational(const SPxRowId& id) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->number(id);
   }



   /// returns index of the column with identifier \p id
   int SoPlex2::idxRational(const SPxColId& id) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->number(id);
   }



   /// returns index of the row or column with identifier \p id
   int SoPlex2::idxRational(const SPxId& id) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->number(id);
   }



   /// gets row \p i
   void SoPlex2::getRowRational(int i, LPRowRational& lprow) const
   {
      assert(_rationalLP != 0);
      _rationalLP->getRow(i, lprow);
   }



   /// gets row with identifier \p id
   void SoPlex2::getRowRational(const SPxRowId& id, LPRowRational& lprow) const
   {
      assert(_rationalLP != 0);
      _rationalLP->getRow(id, lprow);
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



   /// returns vector of row with identifier \p id
   const SVectorRational& SoPlex2::rowVectorRational(const SPxRowId& id) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->rowVector(id);
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



   /// returns right-hand side of row with identifier \p id
   Rational SoPlex2::rhsRational(const SPxRowId& id) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->rhs(id);
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



   /// returns left-hand side of row with identifier \p id
   Rational SoPlex2::lhsRational(const SPxRowId& id) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->lhs(id);
   }



   /// returns inequality type of row \p i
   LPRowRational::Type SoPlex2::rowTypeRational(int i) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->rowType(i);
   }



   /// returns inequality type of row with identifier \p id
   LPRowRational::Type SoPlex2::rowTypeRational(const SPxRowId& id) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->rowType(id);
   }



   /// gets column \p i
   void SoPlex2::getColRational(int i, LPColRational& lpcol) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->getCol(i, lpcol);
   }



   /// gets column with identifier \p id.
   void SoPlex2::getColRational(const SPxColId& id, LPColRational& lpcol) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->getCol(id, lpcol);
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



   /// returns vector of column with identifier \p id
   const SVectorRational& SoPlex2::colVectorRational(const SPxColId& id) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->colVector(id);
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



   /// returns upper bound of column with identifier \p id
   Rational SoPlex2::upperRational(const SPxColId& id) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->upper(id);
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



   /// returns lower bound of column with identifier \p id
   Rational SoPlex2::lowerRational(const SPxColId& id) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->lower(id);
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



   /// returns objective value of column with identifier \p id
   Rational SoPlex2::objRational(const SPxColId& id) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->obj(id);
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



   /// returns objective value of column with identifier \p id after transformation to a maximization problem; since
   /// this is how it is stored internally, this is generally faster
   Rational SoPlex2::maxObjRational(const SPxColId& id) const
   {
      assert(_rationalLP != 0);
      return _rationalLP->maxObj(id);
   }



   /// adds a single row
   void SoPlex2::addRowReal(const LPRowReal& lprow)
   {
      SPxRowId id;
      SoPlex2::addRowReal(id, lprow);
   }



   /// adds a single row and gets its \p id
   void SoPlex2::addRowReal(SPxRowId& id, const LPRowReal& lprow)
   {
      assert(_realLP != 0);
      _realLP->addRow(id, lprow);

      if( _hasBasisReal && !_isRealLPLoaded )
         _basisStatusRowsReal.append(SPxSolver::BASIC);

      _invalidateSolutionReal();
   }



   /// adds multiple rows
   void SoPlex2::addRowsReal(const LPRowSetReal& lprowset)
   {
      assert(_realLP != 0);
      _realLP->addRows(lprowset);

      if( _hasBasisReal && !_isRealLPLoaded )
         _basisStatusRowsReal.append(lprowset.num(), SPxSolver::BASIC);

      _invalidateSolutionReal();
   }



   /// adds multiple rows and gets an array of their \p id 's
   void SoPlex2::addRowsReal(SPxRowId id[], const LPRowSetReal& lprowset)
   {
      assert(_realLP != 0);
      _realLP->addRows(id, lprowset);

      if( _hasBasisReal && !_isRealLPLoaded )
         _basisStatusRowsReal.append(lprowset.num(), SPxSolver::BASIC);

      _invalidateSolutionReal();
   }



   /// adds a single column
   void SoPlex2::addColReal(const LPColReal& lpcol)
   {
      SPxColId id;
      SoPlex2::addColReal(id, lpcol);
   }



   /// adds a single column and gets its \p id
   void SoPlex2::addColReal(SPxColId& id, const LPColReal& lpcol)
   {
      assert(_realLP != 0);
      _realLP->addCol(id, lpcol);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         if( lpcol.lower() > -realParam(SoPlex2::INFTY) )
            _basisStatusColsReal.append(SPxSolver::ON_LOWER);
         else if( lpcol.upper() < realParam(SoPlex2::INFTY) )
            _basisStatusColsReal.append(SPxSolver::ON_UPPER);
         else
            _basisStatusColsReal.append(SPxSolver::ZERO);
      }

      _invalidateSolutionReal();
   }



   /// adds multiple columns
   void SoPlex2::addColsReal(const LPColSetReal& lpcolset)
   {
      assert(_realLP != 0);
      _realLP->addCols(lpcolset);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         for( int i = 0; i < lpcolset.num(); i++ )
         {
            if( lpcolset.lower(i) > -realParam(SoPlex2::INFTY) )
               _basisStatusColsReal.append(SPxSolver::ON_LOWER);
            else if( lpcolset.upper(i) < realParam(SoPlex2::INFTY) )
               _basisStatusColsReal.append(SPxSolver::ON_UPPER);
            else
               _basisStatusColsReal.append(SPxSolver::ZERO);
         }
      }

      _invalidateSolutionReal();
   }



   /// adds multiple columns and gets an array of their \p id 's
   void SoPlex2::addColsReal(SPxColId id[], const LPColSetReal& lpcolset)
   {
      assert(_realLP != 0);
      _realLP->addCols(id, lpcolset);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         for( int i = 0; i < lpcolset.num(); i++ )
         {
            if( lpcolset.lower(i) > -realParam(SoPlex2::INFTY) )
               _basisStatusColsReal.append(SPxSolver::ON_LOWER);
            else if( lpcolset.upper(i) < realParam(SoPlex2::INFTY) )
               _basisStatusColsReal.append(SPxSolver::ON_UPPER);
            else
               _basisStatusColsReal.append(SPxSolver::ZERO);
         }
      }

      _invalidateSolutionReal();
   }



   /// replaces row \p i with \p lprow
   void SoPlex2::changeRowReal(int i, const LPRowReal& lprow)
   {
      assert(_realLP != 0);
      _realLP->changeRow(i, lprow);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         if( _basisStatusRowsReal[i] != SPxSolver::BASIC )
            _hasBasisReal = false;
         else if( _basisStatusRowsReal[i] == SPxSolver::ON_LOWER && lprow.lhs() <= -realParam(SoPlex2::INFTY) )
            _basisStatusRowsReal[i] = (lprow.rhs() < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusRowsReal[i] == SPxSolver::ON_UPPER && lprow.rhs() >= realParam(SoPlex2::INFTY) )
            _basisStatusRowsReal[i] = (lprow.lhs() > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }

      _invalidateSolutionReal();
   }



   /// replaces row with identifier \p id with \p lprow
   void SoPlex2::changeRowReal(SPxRowId id, const LPRowReal& lprow)
   {
      SoPlex2::changeRowReal(idxReal(id), lprow);
   }



   /// changes left-hand side vector for constraints to \p lhs
   void SoPlex2::changeLhsReal(const VectorReal& lhs)
   {
      assert(_realLP != 0);
      _realLP->changeLhs(lhs);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         for( int i = numRowsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusRowsReal[i] == SPxSolver::ON_LOWER && lhs[i] <= -realParam(SoPlex2::INFTY) )
               _basisStatusRowsReal[i] = (rhsReal(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         }
      }

      _invalidateSolutionReal();
   }



   /// changes left-hand side of row \p i to \p lhs
   void SoPlex2::changeLhsReal(int i, Real lhs)
   {
      assert(_realLP != 0);
      _realLP->changeLhs(i, lhs);

      if( _hasBasisReal && !_isRealLPLoaded && _basisStatusRowsReal[i] == SPxSolver::ON_LOWER && lhs <= -realParam(SoPlex2::INFTY) )
         _basisStatusRowsReal[i] = (rhsReal(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;

      _invalidateSolutionReal();
   }



   /// changes left-hand side of row with identifier \p id to \p lhs
   void SoPlex2::changeLhsReal(SPxRowId id, Real lhs)
   {
      SoPlex2::changeLhsReal(idxReal(id), lhs);
   }



   /// changes right-hand side vector to \p rhs
   void SoPlex2::changeRhsReal(const VectorReal& rhs)
   {
      assert(_realLP != 0);
      _realLP->changeRhs(rhs);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         for( int i = numRowsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusRowsReal[i] == SPxSolver::ON_UPPER && rhs[i] >= realParam(SoPlex2::INFTY) )
               _basisStatusRowsReal[i] = (lhsReal(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }

      _invalidateSolutionReal();
   }



   /// changes right-hand side of row \p i to \p rhs
   void SoPlex2::changeRhsReal(int i, Real rhs)
   {
      assert(_realLP != 0);
      _realLP->changeRhs(i, rhs);

      if( _hasBasisReal && !_isRealLPLoaded && _basisStatusRowsReal[i] == SPxSolver::ON_UPPER && rhs >= realParam(SoPlex2::INFTY) )
         _basisStatusRowsReal[i] = (lhsReal(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;

      _invalidateSolutionReal();
   }



   /// changes right-hand of row with identifier \p id to \p rhs
   void SoPlex2::changeRhsReal(SPxRowId id, Real rhs)
   {
      SoPlex2::changeRhsReal(idxReal(id), rhs);
   }



   /// changes left- and right-hand side vectors
   void SoPlex2::changeRangeReal(const VectorReal& lhs, const VectorReal& rhs)
   {
      assert(_realLP != 0);
      _realLP->changeRange(lhs, rhs);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         for( int i = numRowsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusRowsReal[i] == SPxSolver::ON_LOWER && lhs[i] <= -realParam(SoPlex2::INFTY) )
               _basisStatusRowsReal[i] = (rhs[i] < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
            else if( _basisStatusRowsReal[i] == SPxSolver::ON_UPPER && rhs[i] >= realParam(SoPlex2::INFTY) )
               _basisStatusRowsReal[i] = (lhs[i] > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }

      _invalidateSolutionReal();
   }



   /// changes left- and right-hand side of row \p i
   void SoPlex2::changeRangeReal(int i, Real lhs, Real rhs)
   {
      assert(_realLP != 0);
      _realLP->changeRange(i, lhs, rhs);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         if( _basisStatusRowsReal[i] == SPxSolver::ON_LOWER && lhs <= -realParam(SoPlex2::INFTY) )
            _basisStatusRowsReal[i] = (rhs < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusRowsReal[i] == SPxSolver::ON_UPPER && rhs >= realParam(SoPlex2::INFTY) )
            _basisStatusRowsReal[i] = (lhs > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }

      _invalidateSolutionReal();
   }



   /// changes left- and right-hand side of row with identifier \p id
   void SoPlex2::changeRangeReal(SPxRowId id, Real lhs, Real rhs)
   {
      SoPlex2::changeRangeReal(idxReal(id), lhs, rhs);
   }



   /// replaces column \p i with \p lpcol
   void SoPlex2::changeColReal(int i, const LPColReal& lpcol)
   {
      assert(_realLP != 0);
      _realLP->changeCol(i, lpcol);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         if( _basisStatusColsReal[i] == SPxSolver::BASIC )
            _hasBasisReal = false;
         else if( _basisStatusColsReal[i] == SPxSolver::ON_LOWER && lpcol.lower() <= -realParam(SoPlex2::INFTY) )
            _basisStatusColsReal[i] = (lpcol.upper() < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusColsReal[i] == SPxSolver::ON_UPPER && lpcol.upper() >= realParam(SoPlex2::INFTY) )
            _basisStatusColsReal[i] = (lpcol.lower() > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }

      _invalidateSolutionReal();
   }



   /// replaces column with identifier \p id with \p lpcol
   void SoPlex2::changeColReal(SPxColId id, const LPColReal& lpcol)
   {
      SoPlex2::changeColReal(idxReal(id), lpcol);
   }



   /// changes vector of lower bounds to \p lower
   void SoPlex2::changeLowerReal(const VectorReal& lower)
   {
      assert(_realLP != 0);
      _realLP->changeLower(lower);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         for( int i = numColsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusColsReal[i] == SPxSolver::ON_LOWER && lower[i] <= -realParam(SoPlex2::INFTY) )
               _basisStatusColsReal[i] = (upperReal(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         }
      }

      _invalidateSolutionReal();
   }



   /// changes lower bound of column i to \p lower
   void SoPlex2::changeLowerReal(int i, Real lower)
   {
      assert(_realLP != 0);
      _realLP->changeLower(i, lower);

      if( _hasBasisReal && !_isRealLPLoaded && _basisStatusColsReal[i] == SPxSolver::ON_LOWER && lower <= -realParam(SoPlex2::INFTY) )
         _basisStatusColsReal[i] = (upperReal(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;

      _invalidateSolutionReal();
   }



   /// changes lower bound of column with identifier \p id to \p lower
   void SoPlex2::changeLowerReal(SPxColId id, Real lower)
   {
      SoPlex2::changeLowerReal(idxReal(id), lower);
   }



   /// changes vector of upper bounds to \p upper
   void SoPlex2::changeUpperReal(const VectorReal& upper)
   {
      assert(_realLP != 0);
      _realLP->changeUpper(upper);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         for( int i = numColsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusColsReal[i] == SPxSolver::ON_UPPER && upper[i] >= realParam(SoPlex2::INFTY) )
               _basisStatusColsReal[i] = (lowerReal(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }

      _invalidateSolutionReal();
   }



   /// changes \p i 'th upper bound to \p upper
   void SoPlex2::changeUpperReal(int i, Real upper)
   {
      assert(_realLP != 0);
      _realLP->changeUpper(i, upper);

      if( _hasBasisReal && !_isRealLPLoaded &&  _basisStatusColsReal[i] == SPxSolver::ON_UPPER && upper >= realParam(SoPlex2::INFTY) )
         _basisStatusColsReal[i] = (lowerReal(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;

      _invalidateSolutionReal();
   }



   /// changes upper bound of column with identifier \p id to \p upper
   void SoPlex2::changeUpperReal(SPxColId id, Real upper)
   {
      SoPlex2::changeUpperReal(idxReal(id), upper);
   }



   /// changes vectors of column bounds to \p lower and \p upper
   void SoPlex2::changeBoundsReal(const VectorReal& lower, const VectorReal& upper)
   {
      assert(_realLP != 0);
      _realLP->changeBounds(lower, upper);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         for( int i = numColsReal() - 1; i >= 0; i-- )
         {
            if( _basisStatusColsReal[i] == SPxSolver::ON_LOWER && lower[i] <= -realParam(SoPlex2::INFTY) )
               _basisStatusColsReal[i] = (upper[i] < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
            else if( _basisStatusColsReal[i] == SPxSolver::ON_UPPER && upper[i] >= realParam(SoPlex2::INFTY) )
               _basisStatusColsReal[i] = (lower[i] > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }

      _invalidateSolutionReal();
   }



   /// changes bounds of column \p i to \p lower and \p upper
   void SoPlex2::changeBoundsReal(int i, Real lower, Real upper)
   {
      assert(_realLP != 0);
      _realLP->changeBounds(i, lower, upper);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         if( _basisStatusColsReal[i] == SPxSolver::ON_LOWER && lower <= -realParam(SoPlex2::INFTY) )
            _basisStatusColsReal[i] = (upper < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusColsReal[i] == SPxSolver::ON_UPPER && upper >= realParam(SoPlex2::INFTY) )
            _basisStatusColsReal[i] = (lower > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }

      _invalidateSolutionReal();
   }



   /// changes bounds of column with identifier \p id to \p lower and \p upper
   void SoPlex2::changeBoundsReal(SPxColId id, Real lower, Real upper)
   {
      SoPlex2::changeBoundsReal(idxReal(id), lower, upper);
   }



   /// changes objective function vector to \p obj
   void SoPlex2::changeObjReal(const VectorReal& obj)
   {
      assert(_realLP != 0);
      _realLP->changeObj(obj);

      _invalidateSolutionReal();
   }



   /// changes objective coefficient of column i to \p obj
   void SoPlex2::changeObjReal(int i, Real obj)
   {
      assert(_realLP != 0);
      _realLP->changeObj(i, obj);

      _invalidateSolutionReal();
   }



   /// changes objective coefficient of column with identifier \p id to \p obj
   void SoPlex2::changeObjReal(SPxColId id, Real obj)
   {
      assert(_realLP != 0);
      _realLP->changeObj(id, obj);

      _invalidateSolutionReal();
   }



   /// changes matrix entry in row \p i and column \p j to \p val
   void SoPlex2::changeElementReal(int i, int j, Real val)
   {
      assert(_realLP != 0);
      _realLP->changeElement(i, j, val);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         if( _basisStatusRowsReal[i] != SPxSolver::BASIC && _basisStatusColsReal[i] == SPxSolver::BASIC )
            _hasBasisReal = false;
      }

      _invalidateSolutionReal();
   }



   /// changes matrix entry identified by (\p rowid, \p colid) to \p val
   void SoPlex2::changeElementReal(SPxRowId rowid, SPxColId colid, Real val)
   {
      SoPlex2::changeElementReal(idxReal(rowid), idxReal(colid), val);
   }



   /// removes row \p i
   void SoPlex2::removeRowReal(int i)
   {
      assert(_realLP != 0);
      _realLP->removeRow(i);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         if( _basisStatusRowsReal[i] != SPxSolver::BASIC )
            _hasBasisReal = false;
         else
         {
            _basisStatusRowsReal[i] = _basisStatusRowsReal[_basisStatusRowsReal.size() - 1];
            _basisStatusRowsReal.removeLast();
         }
      }

      _invalidateSolutionReal();
   }



   /// removes row with identifier \p id
   void SoPlex2::removeRowReal(SPxRowId id)
   {
      SoPlex2::removeRowReal(idxReal(id));
   }



   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsReal()
   void SoPlex2::removeRowsReal(int perm[])
   {
      assert(_realLP != 0);
      _realLP->removeRows(perm);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         for( int i = numRowsReal() - 1; i >= 0 && _hasBasisReal; i-- )
         {
            if( perm[i] < 0 && _basisStatusRowsReal[i] != SPxSolver::BASIC )
               _hasBasisReal = false;
            else if( perm[i] >= 0 && perm[i] != i )
            {
               assert(perm[i] < numRowsReal());
               assert(perm[perm[i]] < 0);

               _basisStatusRowsReal[perm[i]] = _basisStatusRowsReal[i];
            }
         }

         if( _hasBasisReal )
            _basisStatusRowsReal.reSize(numRowsReal());
      }

      _invalidateSolutionReal();
   }



   /// remove all rows with identifier in array \p id of size \p n; an array \p perm of size #numRowsReal() may be
   /// passed as buffer memory
   void SoPlex2::removeRowsReal(SPxRowId id[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numRowsReal());
         _idToPerm((SPxId*)id, n, p.get_ptr(), numRowsReal());
         SoPlex2::removeRowsReal(p.get_ptr());
      }
      else
      {
         _idToPerm((SPxId*)id, n, perm, numRowsReal());
         SoPlex2::removeRowsReal(perm);
      }
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
      _realLP->removeCol(i);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         if( _basisStatusColsReal[i] == SPxSolver::BASIC )
            _hasBasisReal = false;
         else
         {
            _basisStatusColsReal[i] = _basisStatusColsReal[_basisStatusColsReal.size() - 1];
            _basisStatusColsReal.removeLast();
         }
      }

      _invalidateSolutionReal();
   }



   /// removes column with identifier \p id
   void SoPlex2::removeColReal(SPxColId id)
   {
      SoPlex2::removeColReal(idxReal(id));
   }



   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsReal()
   void SoPlex2::removeColsReal(int perm[])
   {
      assert(_realLP != 0);
      _realLP->removeCols(perm);

      if( _hasBasisReal && !_isRealLPLoaded )
      {
         for( int i = numColsReal() - 1; i >= 0 && _hasBasisReal; i-- )
         {
            if( perm[i] < 0 && _basisStatusColsReal[i] == SPxSolver::BASIC )
               _hasBasisReal = false;
            else if( perm[i] >= 0 && perm[i] != i )
            {
               assert(perm[i] < numColsReal());
               assert(perm[perm[i]] < 0);

               _basisStatusColsReal[perm[i]] = _basisStatusColsReal[i];
            }
         }

         if( _hasBasisReal )
            _basisStatusColsReal.reSize(numColsReal());
      }

      _invalidateSolutionReal();
   }



   /// remove all columns with identifier in array \p id of size \p n; an array \p perm of size #numColsReal() may be
   /// passed as buffer memory
   void SoPlex2::removeColsReal(SPxColId id[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numColsReal());
         _idToPerm((SPxId*)id, n, p.get_ptr(), numColsReal());
         SoPlex2::removeColsReal(p.get_ptr());
      }
      else
      {
         _idToPerm((SPxId*)id, n, perm, numColsReal());
         SoPlex2::removeColsReal(perm);
      }
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
      _hasBasisReal = false;

      _invalidateSolutionReal();
   }



   /// synchronizes real LP with rational LP
   void SoPlex2::syncRealLP()
   {
      assert(_isConsistent());

      // start timing
      _statistics->syncTime.start();

      // copy LP
      if( _isRealLPLoaded )
         _solver.loadLP((SPxLPReal)(*_rationalLP));
      else
         *_realLP = *_rationalLP;

      // load basis if available
      if( _hasBasisRational )
      {
         assert(_basisStatusRowsRational.size() == numRowsReal());
         assert(_basisStatusColsRational.size() == numColsReal());

         if( _isRealLPLoaded )
            _solver.setBasis(_basisStatusRowsRational.get_ptr(), _basisStatusColsRational.get_ptr());
         else
         {
            _basisStatusRowsReal = _basisStatusRowsRational;
            _basisStatusColsReal = _basisStatusColsRational;
         }

         _hasBasisReal = true;
      }
      else
         _hasBasisReal = false;

      // invalidate solution
      _invalidateSolutionReal();

      // stop timing
      _statistics->syncTime.stop();
   }



   /// adds a single row
   void SoPlex2::addRowRational(const LPRowRational& lprow)
   {
      SPxRowId id;
      SoPlex2::addRowRational(id, lprow);
   }



   /// adds a single row and gets its \p id
   void SoPlex2::addRowRational(SPxRowId& id, const LPRowRational& lprow)
   {
      assert(_rationalLP != 0);
      _rationalLP->addRow(id, lprow);

      if( _hasBasisRational )
         _basisStatusRowsRational.append(SPxSolver::BASIC);

      _invalidateSolutionRational();
   }



   /// adds multiple rows
   void SoPlex2::addRowsRational(const LPRowSetRational& lprowset)
   {
      assert(_rationalLP != 0);
      _rationalLP->addRows(lprowset);

      if( _hasBasisRational )
         _basisStatusRowsRational.append(lprowset.num(), SPxSolver::BASIC);

      _invalidateSolutionRational();
   }



   /// adds multiple rows and gets an array of their \p id 's
   void SoPlex2::addRowsRational(SPxRowId id[], const LPRowSetRational& lprowset)
   {
      assert(_rationalLP != 0);
      _rationalLP->addRows(id, lprowset);

      if( _hasBasisRational )
         _basisStatusRowsRational.append(lprowset.num(), SPxSolver::BASIC);

      _invalidateSolutionRational();
   }



   /// adds a single column
   void SoPlex2::addColRational(const LPColRational& lpcol)
   {
      SPxColId id;
      SoPlex2::addColRational(id, lpcol);
   }



   /// adds a single column and gets its \p id
   void SoPlex2::addColRational(SPxColId& id, const LPColRational& lpcol)
   {
      assert(_rationalLP != 0);
      _rationalLP->addCol(id, lpcol);

      if( _hasBasisRational )
      {
         if( lpcol.lower() > -realParam(SoPlex2::INFTY) )
            _basisStatusColsRational.append(SPxSolver::ON_LOWER);
         else if( lpcol.upper() < realParam(SoPlex2::INFTY) )
            _basisStatusColsRational.append(SPxSolver::ON_UPPER);
         else
            _basisStatusColsRational.append(SPxSolver::ZERO);
      }

      _invalidateSolutionRational();
   }



   /// adds multiple columns
   void SoPlex2::addColsRational(const LPColSetRational& lpcolset)
   {
      assert(_rationalLP != 0);
      _rationalLP->addCols(lpcolset);

      if( _hasBasisRational )
      {
         for( int i = 0; i < lpcolset.num(); i++ )
         {
            if( lpcolset.lower(i) > -realParam(SoPlex2::INFTY) )
               _basisStatusColsRational.append(SPxSolver::ON_LOWER);
            else if( lpcolset.upper(i) < realParam(SoPlex2::INFTY) )
               _basisStatusColsRational.append(SPxSolver::ON_UPPER);
            else
               _basisStatusColsRational.append(SPxSolver::ZERO);
         }
      }

      _invalidateSolutionRational();
   }



   /// adds multiple columns and gets an array of their \p id 's
   void SoPlex2::addColsRational(SPxColId id[], const LPColSetRational& lpcolset)
   {
      assert(_rationalLP != 0);
      _rationalLP->addCols(id, lpcolset);

      if( _hasBasisRational )
      {
         for( int i = 0; i < lpcolset.num(); i++ )
         {
            if( lpcolset.lower(i) > -realParam(SoPlex2::INFTY) )
               _basisStatusColsRational.append(SPxSolver::ON_LOWER);
            else if( lpcolset.upper(i) < realParam(SoPlex2::INFTY) )
               _basisStatusColsRational.append(SPxSolver::ON_UPPER);
            else
               _basisStatusColsRational.append(SPxSolver::ZERO);
         }
      }

      _invalidateSolutionRational();
   }



   /// replaces row \p i with \p lprow
   void SoPlex2::changeRowRational(int i, const LPRowRational& lprow)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeRow(i, lprow);

      if( _hasBasisRational )
      {
         if( _basisStatusRowsRational[i] != SPxSolver::BASIC )
            _hasBasisRational = false;
         else if( _basisStatusRowsRational[i] == SPxSolver::ON_LOWER && lprow.lhs() <= -realParam(SoPlex2::INFTY) )
            _basisStatusRowsRational[i] = (lprow.rhs() < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusRowsRational[i] == SPxSolver::ON_UPPER && lprow.rhs() >= realParam(SoPlex2::INFTY) )
            _basisStatusRowsRational[i] = (lprow.lhs() > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }

      _invalidateSolutionRational();
   }



   /// replaces row with identifier \p id with \p lprow
   void SoPlex2::changeRowRational(SPxRowId id, const LPRowRational& lprow)
   {
      SoPlex2::changeRowRational(idxRational(id), lprow);
   }



   /// changes left-hand side vector for constraints to \p lhs
   void SoPlex2::changeLhsRational(const VectorRational& lhs)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeLhs(lhs);

      if( _hasBasisRational )
      {
         for( int i = numRowsRational() - 1; i >= 0; i-- )
         {
            if( _basisStatusRowsRational[i] == SPxSolver::ON_LOWER && lhs[i] <= -realParam(SoPlex2::INFTY) )
               _basisStatusRowsRational[i] = (rhsRational(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         }
      }

      _invalidateSolutionRational();
   }



   /// changes left-hand side of row \p i to \p lhs
   void SoPlex2::changeLhsRational(int i, Rational lhs)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeLhs(i, lhs);

      if( _hasBasisRational && _basisStatusRowsRational[i] == SPxSolver::ON_LOWER && lhs <= -realParam(SoPlex2::INFTY) )
         _basisStatusRowsRational[i] = (rhsRational(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;

      _invalidateSolutionRational();
   }



   /// changes left-hand side of row with identifier \p id to \p lhs
   void SoPlex2::changeLhsRational(SPxRowId id, Rational lhs)
   {
      SoPlex2::changeLhsRational(idxRational(id), lhs);
   }



   /// changes right-hand side vector to \p rhs
   void SoPlex2::changeRhsRational(const VectorRational& rhs)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeRhs(rhs);

      if( _hasBasisRational )
      {
         for( int i = numRowsRational() - 1; i >= 0; i-- )
         {
            if( _basisStatusRowsRational[i] == SPxSolver::ON_UPPER && rhs[i] >= realParam(SoPlex2::INFTY) )
               _basisStatusRowsRational[i] = (lhsRational(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }

      _invalidateSolutionRational();
   }



   /// changes right-hand side of row \p i to \p rhs
   void SoPlex2::changeRhsRational(int i, Rational rhs)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeRhs(i, rhs);

      if( _hasBasisRational && _basisStatusRowsRational[i] == SPxSolver::ON_UPPER && rhs >= realParam(SoPlex2::INFTY) )
         _basisStatusRowsRational[i] = (lhsRational(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;

      _invalidateSolutionRational();
   }



   /// changes right-hand of row with identifier \p id to \p rhs
   void SoPlex2::changeRhsRational(SPxRowId id, Rational rhs)
   {
      SoPlex2::changeRhsRational(idxRational(id), rhs);
   }



   /// changes left- and right-hand side vectors
   void SoPlex2::changeRangeRational(const VectorRational& lhs, const VectorRational& rhs)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeRange(lhs, rhs);

      if( _hasBasisRational )
      {
         for( int i = numRowsRational() - 1; i >= 0; i-- )
         {
            if( _basisStatusRowsRational[i] == SPxSolver::ON_LOWER && lhs[i] <= -realParam(SoPlex2::INFTY) )
               _basisStatusRowsRational[i] = (rhs[i] < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
            else if( _basisStatusRowsRational[i] == SPxSolver::ON_UPPER && rhs[i] >= realParam(SoPlex2::INFTY) )
               _basisStatusRowsRational[i] = (lhs[i] > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }

      _invalidateSolutionRational();
   }



   /// changes left- and right-hand side of row \p i
   void SoPlex2::changeRangeRational(int i, Rational lhs, Rational rhs)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeRange(i, lhs, rhs);

      if( _hasBasisRational )
      {
         if( _basisStatusRowsRational[i] == SPxSolver::ON_LOWER && lhs <= -realParam(SoPlex2::INFTY) )
            _basisStatusRowsRational[i] = (rhs < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusRowsRational[i] == SPxSolver::ON_UPPER && rhs >= realParam(SoPlex2::INFTY) )
            _basisStatusRowsRational[i] = (lhs > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }

      _invalidateSolutionRational();
   }



   /// changes left- and right-hand side of row with identifier \p id
   void SoPlex2::changeRangeRational(SPxRowId id, Rational lhs, Rational rhs)
   {
      SoPlex2::changeRangeRational(idxRational(id), lhs, rhs);
   }



   /// replaces column \p i with \p lpcol
   void SoPlex2::changeColRational(int i, const LPColRational& lpcol)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeCol(i, lpcol);

      if( _hasBasisRational )
      {
         if( _basisStatusColsRational[i] == SPxSolver::BASIC )
            _hasBasisRational = false;
         else if( _basisStatusColsRational[i] == SPxSolver::ON_LOWER && lpcol.lower() <= -realParam(SoPlex2::INFTY) )
            _basisStatusColsRational[i] = (lpcol.upper() < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusColsRational[i] == SPxSolver::ON_UPPER && lpcol.upper() >= realParam(SoPlex2::INFTY) )
            _basisStatusColsRational[i] = (lpcol.lower() > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }

      _invalidateSolutionRational();
   }



   /// replaces column with identifier \p id with \p lpcol
   void SoPlex2::changeColRational(SPxColId id, const LPColRational& lpcol)
   {
      SoPlex2::changeColRational(idxRational(id), lpcol);
   }



   /// changes vector of lower bounds to \p lower
   void SoPlex2::changeLowerRational(const VectorRational& lower)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeLower(lower);

      if( _hasBasisRational )
      {
         for( int i = numColsRational() - 1; i >= 0; i-- )
         {
            if( _basisStatusColsRational[i] == SPxSolver::ON_LOWER && lower[i] <= -realParam(SoPlex2::INFTY) )
               _basisStatusColsRational[i] = (upperRational(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         }
      }

      _invalidateSolutionRational();
   }



   /// changes lower bound of column i to \p lower
   void SoPlex2::changeLowerRational(int i, Rational lower)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeLower(i, lower);

      if( _hasBasisRational && _basisStatusColsRational[i] == SPxSolver::ON_LOWER && lower <= -realParam(SoPlex2::INFTY) )
         _basisStatusColsRational[i] = (upperRational(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;

      _invalidateSolutionRational();
   }



   /// changes lower bound of column with identifier \p id to \p lower
   void SoPlex2::changeLowerRational(SPxColId id, Rational lower)
   {
      SoPlex2::changeLowerRational(idxRational(id), lower);
   }



   /// changes vector of upper bounds to \p upper
   void SoPlex2::changeUpperRational(const VectorRational& upper)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeUpper(upper);

      if( _hasBasisRational )
      {
         for( int i = numColsRational() - 1; i >= 0; i-- )
         {
            if( _basisStatusColsRational[i] == SPxSolver::ON_UPPER && upper[i] >= realParam(SoPlex2::INFTY) )
               _basisStatusColsRational[i] = (lowerRational(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }

      _invalidateSolutionRational();
   }



   /// changes \p i 'th upper bound to \p upper
   void SoPlex2::changeUpperRational(int i, Rational upper)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeUpper(i, upper);

      if( _hasBasisRational &&  _basisStatusColsRational[i] == SPxSolver::ON_UPPER && upper >= realParam(SoPlex2::INFTY) )
         _basisStatusColsRational[i] = (lowerRational(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;

      _invalidateSolutionRational();
   }



   /// changes upper bound of column with identifier \p id to \p upper
   void SoPlex2::changeUpperRational(SPxColId id, Rational upper)
   {
      SoPlex2::changeUpperRational(idxRational(id), upper);
   }



   /// changes vectors of column bounds to \p lower and \p upper
   void SoPlex2::changeBoundsRational(const VectorRational& lower, const VectorRational& upper)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeBounds(lower, upper);

      if( _hasBasisRational )
      {
         for( int i = numColsRational() - 1; i >= 0; i-- )
         {
            if( _basisStatusColsRational[i] == SPxSolver::ON_LOWER && lower[i] <= -realParam(SoPlex2::INFTY) )
               _basisStatusColsRational[i] = (upper[i] < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
            else if( _basisStatusColsRational[i] == SPxSolver::ON_UPPER && upper[i] >= realParam(SoPlex2::INFTY) )
               _basisStatusColsRational[i] = (lower[i] > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
         }
      }

      _invalidateSolutionRational();
   }



   /// changes bounds of column \p i to \p lower and \p upper
   void SoPlex2::changeBoundsRational(int i, Rational lower, Rational upper)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeBounds(i, lower, upper);

      if( _hasBasisRational )
      {
         if( _basisStatusColsRational[i] == SPxSolver::ON_LOWER && lower <= -realParam(SoPlex2::INFTY) )
            _basisStatusColsRational[i] = (upper < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
         else if( _basisStatusColsRational[i] == SPxSolver::ON_UPPER && upper >= realParam(SoPlex2::INFTY) )
            _basisStatusColsRational[i] = (lower > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
      }

      _invalidateSolutionRational();
   }



   /// changes bounds of column with identifier \p id to \p lower and \p upper
   void SoPlex2::changeBoundsRational(SPxColId id, Rational lower, Rational upper)
   {
      SoPlex2::changeBoundsRational(idxRational(id), lower, upper);
   }



   /// changes objective function vector to \p obj
   void SoPlex2::changeObjRational(const VectorRational& obj)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeObj(obj);

      _invalidateSolutionRational();
   }



   /// changes objective coefficient of column i to \p obj
   void SoPlex2::changeObjRational(int i, Rational obj)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeObj(i, obj);

      _invalidateSolutionRational();
   }



   /// changes objective coefficient of column with identifier \p id to \p obj
   void SoPlex2::changeObjRational(SPxColId id, Rational obj)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeObj(id, obj);

      _invalidateSolutionRational();
   }



   /// changes matrix entry in row \p i and column \p j to \p val
   void SoPlex2::changeElementRational(int i, int j, Rational val)
   {
      assert(_rationalLP != 0);
      _rationalLP->changeElement(i, j, val);

      if( _hasBasisRational )
      {
         if( _basisStatusRowsRational[i] != SPxSolver::BASIC && _basisStatusColsRational[i] == SPxSolver::BASIC )
            _hasBasisRational = false;
      }

      _invalidateSolutionRational();
   }



   /// changes matrix entry identified by (\p rowid, \p colid) to \p val
   void SoPlex2::changeElementRational(SPxRowId rowid, SPxColId colid, Rational val)
   {
      SoPlex2::changeElementRational(idxRational(rowid), idxRational(colid), val);
   }



   /// removes row \p i
   void SoPlex2::removeRowRational(int i)
   {
      assert(_rationalLP != 0);
      _rationalLP->removeRow(i);

      if( _hasBasisRational )
      {
         if( _basisStatusRowsRational[i] != SPxSolver::BASIC )
            _hasBasisRational = false;
         else
         {
            _basisStatusRowsRational[i] = _basisStatusRowsRational[_basisStatusRowsRational.size() - 1];
            _basisStatusRowsRational.removeLast();
         }
      }

      _invalidateSolutionRational();
   }



   /// removes row with identifier \p id
   void SoPlex2::removeRowRational(SPxRowId id)
   {
      SoPlex2::removeRowRational(idxRational(id));
   }



   /// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the new
   /// index where row \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numRowsRational()
   void SoPlex2::removeRowsRational(int perm[])
   {
      assert(_rationalLP != 0);
      _rationalLP->removeRows(perm);

      if( _hasBasisRational )
      {
         for( int i = numRowsRational() - 1; i >= 0 && _hasBasisRational; i-- )
         {
            if( perm[i] < 0 && _basisStatusRowsRational[i] != SPxSolver::BASIC )
               _hasBasisRational = false;
            else if( perm[i] >= 0 && perm[i] != i )
            {
               assert(perm[i] < numRowsRational());
               assert(perm[perm[i]] < 0);

               _basisStatusRowsRational[perm[i]] = _basisStatusRowsRational[i];
            }
         }

         if( _hasBasisRational )
            _basisStatusRowsRational.reSize(numRowsRational());
      }

      _invalidateSolutionRational();
   }



   /// remove all rows with identifier in array \p id of size \p n; an array \p perm of size #numRowsRational() may be
   /// passed as buffer memory
   void SoPlex2::removeRowsRational(SPxRowId id[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numRowsRational());
         _idToPerm((SPxId*)id, n, p.get_ptr(), numRowsRational());
         SoPlex2::removeRowsRational(p.get_ptr());
      }
      else
      {
         _idToPerm((SPxId*)id, n, perm, numRowsRational());
         SoPlex2::removeRowsRational(perm);
      }
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
      _rationalLP->removeCol(i);

      if( _hasBasisRational )
      {
         if( _basisStatusColsRational[i] == SPxSolver::BASIC )
            _hasBasisRational = false;
         else
         {
            _basisStatusColsRational[i] = _basisStatusColsRational[_basisStatusColsRational.size() - 1];
            _basisStatusColsRational.removeLast();
         }
      }

      _invalidateSolutionRational();
   }



   /// removes column with identifier \p id
   void SoPlex2::removeColRational(SPxColId id)
   {
      SoPlex2::removeColRational(idxRational(id));
   }



   /// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
   /// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
   /// #numColsRational()
   void SoPlex2::removeColsRational(int perm[])
   {
      assert(_rationalLP != 0);
      _rationalLP->removeCols(perm);

      if( _hasBasisRational )
      {
         for( int i = numColsRational() - 1; i >= 0 && _hasBasisRational; i-- )
         {
            if( perm[i] < 0 && _basisStatusColsRational[i] == SPxSolver::BASIC )
               _hasBasisRational = false;
            else if( perm[i] >= 0 && perm[i] != i )
            {
               assert(perm[i] < numColsRational());
               assert(perm[perm[i]] < 0);

               _basisStatusColsRational[perm[i]] = _basisStatusColsRational[i];
            }
         }

         if( _hasBasisRational )
            _basisStatusColsRational.reSize(numColsRational());
      }

      _invalidateSolutionRational();
   }



   /// remove all columns with identifier in array \p id of size \p n; an array \p perm of size #numColsRational() may
   /// be passed as buffer memory
   void SoPlex2::removeColsRational(SPxColId id[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numColsRational());
         _idToPerm((SPxId*)id, n, p.get_ptr(), numColsRational());
         SoPlex2::removeColsRational(p.get_ptr());
      }
      else
      {
         _idToPerm((SPxId*)id, n, perm, numColsRational());
         SoPlex2::removeColsRational(perm);
      }
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

      _rationalLP->clear();
      _hasBasisRational = false;

      _invalidateSolutionRational();
   }



   /// synchronizes rational LP with real LP
   void SoPlex2::syncRationalLP()
   {
      assert(_isConsistent());

      // start timing
      _statistics->syncTime.start();

      // copy LP
      *_rationalLP = *_realLP;

      // load basis if available
      if( _hasBasisReal )
      {
         assert(_basisStatusRowsReal.size() == numRowsRational());
         assert(_basisStatusColsReal.size() == numColsRational());

         _basisStatusRowsRational = _basisStatusRowsReal;
         _basisStatusColsRational = _basisStatusColsReal;

         _hasBasisRational = true;
      }
      else
         _hasBasisRational = false;

      // invalidate solution
      _invalidateSolutionRational();

      // stop timing
      _statistics->syncTime.stop();
   }



   /// solves real LP
   SPxSolver::Status SoPlex2::solveReal()
   {
      assert(_isConsistent());

      // clear statistics
      _statistics->clearSolvingData();

      // start timing
      _statistics->solvingTime.start();

      // the solution is no longer valid
      _invalidateSolutionReal();

      // will preprocessing be applied? (only if no basis is available)
      bool applyPreprocessing = (_firstScaler != 0 || _simplifier != 0 || _secondScaler != 0) && !_hasBasisReal;

      if( applyPreprocessing )
      {
         _enableSimplifierAndScalers();
         _solver.setTerminationValue(realParam(SoPlex2::INFTY));
      }
      else
      {
         _disableSimplifierAndScalers();
         ///@todo implement for both objective senses
         _solver.setTerminationValue(intParam(SoPlex2::OBJSENSE) == SoPlex2::OBJSENSE_MINIMIZE
            ? realParam(SoPlex2::OBJLIMIT_UPPER) : realParam(SoPlex2::INFTY));
      }

      if( _isRealLPLoaded )
      {
         assert(_realLP == &_solver);
         assert(!applyPreprocessing || !_hasBasisReal);

         // preprocessing is always applied to the LP in the solver; hence we have to create a copy of the original LP
         // if preprocessing is turned on
         if( applyPreprocessing )
         {
            _realLP = 0;
            spx_alloc(_realLP);
            _realLP = new (_realLP) SPxLPReal(_solver);
            _isRealLPLoaded = false;
         }
         ///@todo maybe this should move closer to the actual solving routine
         else if( _hasBasisReal )
         {
            assert(_solver.basis().status() >= SPxBasis::REGULAR);

            _basisStatusRowsReal.reSize(numRowsReal());
            _basisStatusColsReal.reSize(numColsReal());
            _solver.getBasis(_basisStatusRowsReal.get_ptr(), _basisStatusColsReal.get_ptr());
         }
      }
      else
      {
         assert(_realLP != &_solver);
         assert(!applyPreprocessing || !_hasBasisReal);

         // ensure that the solver has the original problem
         _solver.loadLP(*_realLP);

         // load basis if available
         if( _hasBasisReal )
         {
            assert(_basisStatusRowsReal.size() == numRowsReal());
            assert(_basisStatusColsReal.size() == numColsReal());

            ///@todo this should not fail even if the basis is invalid (wrong dimension or wrong number of basic
            ///      entries); fix either in SPxSolver or in SPxBasis
            _solver.setBasis(_basisStatusRowsReal.get_const_ptr(), _basisStatusColsReal.get_const_ptr());
         }

         // if there is no preprocessing, then the original and the transformed problem are identical and it is more
         // memory-efficient to keep only the problem in the solver
         if( !applyPreprocessing )
         {
            spx_free(_realLP);
            _realLP = &_solver;
            _isRealLPLoaded = true;
         }
      }

      // assert that we have two problems if and only if we apply preprocessing
      assert(_realLP == &_solver || applyPreprocessing);
      assert(_realLP != &_solver || !applyPreprocessing);

      // apply scaling before the simplification
      if( _firstScaler != 0 )
      {
         _firstScaler->scale(_solver);
      }

      // apply problem simplification
      SPxSimplifier::Result simplificationStatus = SPxSimplifier::OKAY;
      if( _simplifier != 0 )
      {
         simplificationStatus = _simplifier->simplify(_solver, realParam(SoPlex2::EPSILON_ZERO), rationalParam(SoPlex2::FEASTOL), rationalParam(SoPlex2::OPTTOL));
      }

      // apply scaling after the simplification
      if( _secondScaler != 0 && simplificationStatus == SPxSimplifier::OKAY )
      {
         _secondScaler->scale(_solver);
      }

      // run the simplex method if problem has not been solved by the simplifier
      if( simplificationStatus == SPxSimplifier::OKAY )
      {
         ///@todo this should be a separate method implementing starter, auto pricing, and recovery mechanisms
         _solver.solve();
      }

      ///@todo move to private helper methods
      // evaluate status flag
      assert(_statusReal == SPxSolver::UNKNOWN);
      if( simplificationStatus == SPxSimplifier::INFEASIBLE )
         _statusReal = SPxSolver::INFEASIBLE;
      else if( simplificationStatus == SPxSimplifier::DUAL_INFEASIBLE )
         _statusReal = SPxSolver::INForUNBD;
      else if( simplificationStatus == SPxSimplifier::UNBOUNDED )
         _statusReal = SPxSolver::UNBOUNDED;
      else if( simplificationStatus == SPxSimplifier::VANISHED )
         _statusReal = SPxSolver::OPTIMAL;
      else if( simplificationStatus == SPxSimplifier::OKAY )
      {
         _statusReal = _solver.status();

         ///@todo move to private helper methods
         // evaluate solution flags
         assert(!_hasPrimalReal);
         assert(_solver.basis().status() != SPxBasis::PRIMAL || statusReal() != SPxSolver::ERROR);
         assert(_solver.basis().status() != SPxBasis::PRIMAL || statusReal() != SPxSolver::NO_RATIOTESTER);
         assert(_solver.basis().status() != SPxBasis::PRIMAL || statusReal() != SPxSolver::NO_PRICER);
         assert(_solver.basis().status() != SPxBasis::PRIMAL || statusReal() != SPxSolver::NO_SOLVER);
         assert(_solver.basis().status() != SPxBasis::PRIMAL || statusReal() != SPxSolver::NOT_INIT);
         assert(_solver.basis().status() != SPxBasis::PRIMAL || statusReal() != SPxSolver::SINGULAR);
         assert(_solver.basis().status() != SPxBasis::PRIMAL || statusReal() != SPxSolver::NO_PROBLEM);
         assert(_solver.basis().status() != SPxBasis::PRIMAL || statusReal() != SPxSolver::UNBOUNDED);
         assert(_solver.basis().status() != SPxBasis::PRIMAL || statusReal() != SPxSolver::INFEASIBLE);
         assert(_solver.basis().status() != SPxBasis::UNBOUNDED || statusReal() == SPxSolver::UNBOUNDED);
         assert(_solver.basis().status() == SPxBasis::UNBOUNDED || statusReal() != SPxSolver::UNBOUNDED);
         _hasPrimalReal = (SPxSolver::OPTIMAL || (_simplifier == 0 &&
               (_solver.basis().status() == SPxBasis::PRIMAL || _solver.basis().status() == SPxBasis::UNBOUNDED) &&
               _solver.shift() < 10.0 * realParam(SoPlex2::EPSILON_ZERO)));

         assert(!_hasPrimalrayReal);
         _hasPrimalrayReal = (statusReal() == SPxSolver::UNBOUNDED && _simplifier == 0);

         assert(!_hasDualReal);
         assert(_solver.basis().status() != SPxBasis::DUAL || statusReal() != SPxSolver::ERROR);
         assert(_solver.basis().status() != SPxBasis::DUAL || statusReal() != SPxSolver::NO_RATIOTESTER);
         assert(_solver.basis().status() != SPxBasis::DUAL || statusReal() != SPxSolver::NO_PRICER);
         assert(_solver.basis().status() != SPxBasis::DUAL || statusReal() != SPxSolver::NO_SOLVER);
         assert(_solver.basis().status() != SPxBasis::DUAL || statusReal() != SPxSolver::NOT_INIT);
         assert(_solver.basis().status() != SPxBasis::DUAL || statusReal() != SPxSolver::SINGULAR);
         assert(_solver.basis().status() != SPxBasis::DUAL || statusReal() != SPxSolver::NO_PROBLEM);
         assert(_solver.basis().status() != SPxBasis::DUAL || statusReal() != SPxSolver::UNBOUNDED);
         assert(_solver.basis().status() != SPxBasis::DUAL || statusReal() != SPxSolver::INFEASIBLE);
         assert(_solver.basis().status() != SPxBasis::INFEASIBLE || statusReal() == SPxSolver::INFEASIBLE);
         assert(_solver.basis().status() == SPxBasis::INFEASIBLE || statusReal() != SPxSolver::INFEASIBLE);
         _hasDualReal = (statusReal() == SPxSolver::OPTIMAL || (_simplifier == 0 &&
               (_solver.basis().status() == SPxBasis::DUAL || _solver.basis().status() == SPxBasis::INFEASIBLE) &&
               _solver.shift() < 10.0 * realParam(SoPlex2::EPSILON_ZERO)));

         assert(!_hasDualfarkasReal);
         _hasDualfarkasReal = (statusReal() == SPxSolver::INFEASIBLE && _simplifier == 0);

         ///@todo move to private helper methods
         // process result
         switch( statusReal() )
         {
         case SPxSolver::OPTIMAL:
            // unsimplify if simplifier is active and LP is solved to optimality; this must be done here and not at solution
            // query, because we want to have the basis for the original problem
            if( _simplifier != 0 )
            {
               assert(!_simplifier->isUnsimplified());
               assert(simplificationStatus == SPxSimplifier::VANISHED || simplificationStatus == SPxSimplifier::OKAY);
               assert(!_hasBasisReal);

               bool vanished = simplificationStatus == SPxSimplifier::VANISHED;

               // get solution vectors for transformed problem
               DVectorReal primal(vanished ? 0 : _solver.nCols());
               DVectorReal slacks(vanished ? 0 : _solver.nRows());
               DVectorReal dual(vanished ? 0 : _solver.nRows());
               DVectorReal redcost(vanished ? 0 : _solver.nCols());

               _basisStatusRowsReal.reSize(numRowsReal());
               _basisStatusColsReal.reSize(numColsReal());
               assert(vanished || _basisStatusRowsReal.size() >= _solver.nRows());
               assert(vanished || _basisStatusColsReal.size() >= _solver.nCols());

               if( !vanished )
               {
                  assert(_solver.status() == SPxSolver::OPTIMAL);

                  _solver.getPrimal(primal);
                  _solver.getSlacks(slacks);
                  _solver.getDual(dual);
                  _solver.getRedCost(redcost);

                  // unscale vectors w.r.t. second scaler
                  if( _secondScaler != 0 )
                  {
                     _secondScaler->unscalePrimal(primal);
                     _secondScaler->unscaleSlacks(slacks);
                     _secondScaler->unscaleDual(dual);
                     _secondScaler->unscaleRedCost(redcost);
                  }

                  // get basis of transformed problem
                  _solver.getBasis(_basisStatusRowsReal.get_ptr(), _basisStatusColsReal.get_ptr());
               }

               ///@todo catch exception
               _simplifier->unsimplify(primal, dual, slacks, redcost, _basisStatusRowsReal.get_ptr(), _basisStatusColsReal.get_ptr());

               // store basis for original problem
               _simplifier->getBasis(_basisStatusRowsReal.get_ptr(), _basisStatusColsReal.get_ptr());
            }
            // if the original problem is not in the solver because of scaling, we also need to store the basis
            else if( !_isRealLPLoaded )
            {
               _basisStatusRowsReal.reSize(numRowsReal());
               _basisStatusColsReal.reSize(numColsReal());
               assert(_basisStatusRowsReal.size() == _solver.nRows());
               assert(_basisStatusColsReal.size() == _solver.nCols());

               _solver.getBasis(_basisStatusRowsReal.get_ptr(), _basisStatusColsReal.get_ptr());
            }

            // in all cases we have a basis for warmstarting
            _hasBasisReal = true;
            break;

         case SPxSolver::ABORT_CYCLING:
         case SPxSolver::ABORT_TIME:
         case SPxSolver::ABORT_ITER:
         case SPxSolver::ABORT_VALUE:
         case SPxSolver::REGULAR:
         case SPxSolver::RUNNING:
         case SPxSolver::UNBOUNDED:
         case SPxSolver::INFEASIBLE:
         case SPxSolver::INForUNBD:
            // store basis if it is regular and the original problem is not in the solver because of scaling
            if( _simplifier == 0 && !_isRealLPLoaded )
            {
               _basisStatusRowsReal.reSize(numRowsReal());
               _basisStatusColsReal.reSize(numColsReal());
               assert(_basisStatusRowsReal.size() == _solver.nRows());
               assert(_basisStatusColsReal.size() == _solver.nCols());

               _solver.getBasis(_basisStatusRowsReal.get_ptr(), _basisStatusColsReal.get_ptr());
               _hasBasisReal = true;
            }
            // non-optimal basis should currently not be unsimplified
            else
               _hasBasisReal = false;
            break;

         case SPxSolver::SINGULAR:
            // if there was a regular starting basis and the original problem is in the solver, load the basis
            if( _hasBasisReal && _isRealLPLoaded )
            {
               assert(_simplifier == 0);
               assert(_basisStatusRowsReal.size() == _solver.nRows());
               assert(_basisStatusColsReal.size() == _solver.nCols());
               _solver.setBasis(_basisStatusRowsReal.get_ptr(), _basisStatusColsReal.get_ptr());
            }
            break;

         default:
            _hasBasisReal = false;
            break;
         }
      }

      // stop timing
      _statistics->solvingTime.stop();

      return statusReal();
   }



   /// returns the current status
   SPxSolver::Status SoPlex2::statusReal() const
   {
      return _statusReal;
   }



   /// returns the current basis status
   SPxBasis::SPxStatus SoPlex2::basisStatusReal() const
   {
      return _solver.basis().status();
   }



   /// returns the objective value if a primal solution is available
   Real SoPlex2::objValueReal() const
   {
      assert(OBJSENSE_MAXIMIZE == 1);
      assert(OBJSENSE_MINIMIZE == -1);

      if( hasPrimalReal() )
      {
         ///@todo remember if computed once
         DVectorReal primal(numColsReal());
         getPrimalReal(primal);
         return (primal * maxObjReal()) * intParam(SoPlex2::OBJSENSE);
      }
      else
         return -realParam(SoPlex2::INFTY) * intParam(SoPlex2::OBJSENSE);
   }



   /// returns the termination value
   Real SoPlex2::terminationValueReal() const
   {
      return _solver.terminationValue();
   }



   /// is a primal feasible solution available?
   bool SoPlex2::hasPrimalReal() const
   {
      return _hasPrimalReal;
   }



   /// gets the primal solution vector if available; returns true on success
   bool SoPlex2::getPrimalReal(VectorReal& vector) const
   {
      if( hasPrimalReal() )
      {
         if( _simplifier != 0 )
         {
            assert(_simplifier->isUnsimplified());
            vector = _simplifier->unsimplifiedPrimal();
         }
         else
         {
            _solver.getPrimal(vector);

            if( _secondScaler != 0 )
               _secondScaler->unscalePrimal(vector);
         }

         if( _firstScaler != 0 )
            _firstScaler->unscalePrimal(vector);

         return true;
      }
      else
         return false;
   }



   /// gets the vector of slack values if available; returns true on success
   bool SoPlex2::getSlacksReal(VectorReal& vector) const
   {
      if( hasPrimalReal() )
      {
         if( _simplifier != 0 )
         {
            assert(_simplifier->isUnsimplified());
            vector = _simplifier->unsimplifiedSlacks();
         }
         else
         {
            _solver.getSlacks(vector);

            if( _secondScaler != 0 )
               _secondScaler->unscaleSlacks(vector);
         }

         if( _firstScaler != 0 )
            _firstScaler->unscaleSlacks(vector);

         return true;
      }
      else
         return false;
   }



   /// is a primal unbounded ray available?
   bool SoPlex2::hasPrimalrayReal() const
   {
      return _hasPrimalrayReal;
   }



   /// gets the primal ray if available; returns true on success
   bool SoPlex2::getPrimalrayReal(VectorReal& vector) const
   {
      if( hasPrimalrayReal() )
      {
         _solver.getPrimalray(vector);

         if( _secondScaler != 0 )
            _secondScaler->unscalePrimal(vector);

         if( _firstScaler != 0 )
            _firstScaler->unscalePrimal(vector);

         return true;
      }
      else
         return false;
   }



   /// is a dual feasible solution available?
   bool SoPlex2::hasDualReal() const
   {
      return _hasDualReal;
   }



   /// gets the dual solution vector if available; returns true on success
   bool SoPlex2::getDualReal(VectorReal& vector) const
   {
      if( hasDualReal() )
      {
         if( _simplifier != 0 )
         {
            assert(_simplifier->isUnsimplified());
            vector = _simplifier->unsimplifiedDual();
         }
         else
         {
            _solver.getDual(vector);

            if( _secondScaler != 0 )
               _secondScaler->unscaleDual(vector);
         }

         if( _firstScaler != 0 )
            _firstScaler->unscaleDual(vector);

         return true;
      }
      else
         return false;
   }



   /// gets the vector of reduced cost values if available; returns true on success
   bool SoPlex2::getRedcostReal(VectorReal& vector) const
   {
      if( hasDualReal() )
      {
         if( _simplifier != 0 )
         {
            assert(_simplifier->isUnsimplified());
            vector = _simplifier->unsimplifiedRedCost();
         }
         else
         {
            _solver.getRedCost(vector);

            if( _secondScaler != 0 )
               _secondScaler->unscaleRedCost(vector);
         }

         if( _firstScaler != 0 )
            _firstScaler->unscaleRedCost(vector);

         return true;
      }
      else
         return false;
   }



   /// is Farkas proof of infeasibility available?
   bool SoPlex2::hasDualfarkasReal() const
   {
      return _hasDualfarkasReal;
   }



   /// gets the Farkas proof if available; returns true on success
   bool SoPlex2::getDualfarkasReal(VectorReal& vector) const
   {
      if( hasDualfarkasReal() )
      {
         _solver.getDualfarkas(vector);

         if( _secondScaler != 0 )
            _secondScaler->unscaleDual(vector);

         if( _firstScaler != 0 )
            _firstScaler->unscaleDual(vector);

         return true;
      }
      else
         return false;
   }



   /// gets violation of bounds by given primal solution
   void SoPlex2::getBoundViolationReal(VectorReal& primal, Real& maxviol, Real& sumviol) const
   {
      assert(primal.dim() >= numColsReal());

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



   /// gets internal violation of bounds by given primal solution
   void SoPlex2::getInternalBoundViolationReal(Real& maxviol, Real& sumviol) const
   {
      maxviol = 0.0;
      sumviol = 0.0;

      DVectorReal primal(_solver.nCols());

      _solver.getPrimal(primal);

      for( int i = _solver.nCols() - 1; i >= 0; i-- )
      {
         Real viol = _solver.lower(i) - primal[i];
         if( viol > 0.0 )
         {
            sumviol += viol;
            if( viol > maxviol )
               maxviol = viol;
         }

         viol = primal[i] - _solver.upper(i);
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



   /// gets internal violation of constraints
   void SoPlex2::getInternalConstraintViolationReal(Real& maxviol, Real& sumviol) const
   {
      _solver.qualConstraintViolation(maxviol, sumviol);
   }



   /// gets violation of slacks
   void SoPlex2::getSlackViolationReal(Real& maxviol, Real& sumviol) const
   {
      _solver.qualSlackViolation(maxviol, sumviol);
   }



   /// gets violation of reduced costs
   void SoPlex2::getRedCostViolationReal(Real& maxviol, Real& sumviol) const
   {
      _solver.qualRedCostViolation(maxviol, sumviol);
   }



   /// synchronizes LPs, clears statistics, and solves rational LP
   SPxSolver::Status SoPlex2::solveRational()
   {
      assert(_isConsistent());
      assert(_statistics != 0);

      // clear statistics
      _statistics->clearSolvingData();

      // start timing
      _statistics->solvingTime.start();

#if 1
      // copy rounded rational LP to real LP
      syncRealLP();

      // call rational solving routine
      _solveRational();
#else
      // copy rounded rational LP to real LP
      syncRealLP();

      // solve floating-point LP
      _statusRational = solveReal();

      // store real as rational solution
      _syncRationalSolution(true, true, true);
#endif

      // stop timing
      _statistics->solvingTime.stop();

      return _statusRational;
   }



   /// returns the current status
   SPxSolver::Status SoPlex2::statusRational() const
   {
      return _statusRational;
   }



   /// returns the objective value if a primal solution is available
   Rational SoPlex2::objValueRational() const
   {
      assert(OBJSENSE_MAXIMIZE == 1);
      assert(OBJSENSE_MINIMIZE == -1);

      if( hasPrimalrayRational() )
      {
         return realParam(SoPlex2::INFTY) * intParam(SoPlex2::OBJSENSE);
      }
      else if( hasPrimalRational() )
      {
         ///@todo remember if computed once
         DVectorRational primal(numColsRational());
         getPrimalRational(primal);
         return (primal * maxObjRational()) * (Rational)intParam(SoPlex2::OBJSENSE);
      }
      else
         return -realParam(SoPlex2::INFTY) * intParam(SoPlex2::OBJSENSE);
   }



   /// is a primal feasible solution available?
   bool SoPlex2::hasPrimalRational() const
   {
      return _solRational.hasPrimal();
   }



   /// gets the primal solution vector if available; returns true on success
   bool SoPlex2::getPrimalRational(VectorRational& vector) const
   {
      return _solRational.getPrimal(vector);
   }



   /// gets the vector of slack values if available; returns true on success
   bool SoPlex2::getSlacksRational(VectorRational& vector) const
   {
      return _solRational.getSlacks(vector);
   }



   /// is a primal unbounded ray available?
   bool SoPlex2::hasPrimalrayRational() const
   {
      return _solRational.hasPrimalray();
   }



   /// gets the primal ray if LP is unbounded; returns true on success
   bool SoPlex2::getPrimalrayRational(VectorRational& vector) const
   {
      return _solRational.getPrimalray(vector);
   }



   /// is a dual feasible solution available?
   bool SoPlex2::hasDualRational() const
   {
      return _solRational.hasDual();
   }



   /// gets the dual solution vector if available; returns true on success
   bool SoPlex2::getDualRational(VectorRational& vector) const
   {
      return _solRational.getDual(vector);
   }



   /// gets the vector of reduced cost values if available; returns true on success
   bool SoPlex2::getRedcostRational(VectorRational& vector) const
   {
      return _solRational.getRedcost(vector);
   }



   /// is Farkas proof of infeasibility available?
   bool SoPlex2::hasDualfarkasRational() const
   {
      return _solRational.hasDualfarkas();
   }

   /// gets the Farkas proof if LP is infeasible; returns true on success
   bool SoPlex2::getDualfarkasRational(VectorRational& vector) const
   {
      return _solRational.getDualfarkas(vector);
   }



   /// gets violation of bounds by given primal solution
   void SoPlex2::getBoundViolationRational(VectorRational& primal, Rational& maxviol, Rational& sumviol) const
   {
      assert(primal.dim() >= numColsRational());

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
   void SoPlex2::getSlackViolationRational(Rational& maxviol, Rational& sumviol) const
   {
      ///@todo implement
      maxviol = 0;
      sumviol = 0;
   }



   /// gets violation of reduced costs
   void SoPlex2::getRedCostViolationRational(Rational& maxviol, Rational& sumviol) const
   {
      ///@todo implement
      maxviol = 0;
      sumviol = 0;
   }



   /// is an advanced starting basis available?
   bool SoPlex2::hasBasisReal() const
   {
      return _hasBasisReal;
   }



   /// returns basis status for a single row
   SPxSolver::VarStatus SoPlex2::basisRowStatusReal(int row) const
   {
      assert(row >= 0);
      assert(row < numRowsReal());

      // if no basis is available, return slack basis
      if( !_hasBasisReal )
         return SPxSolver::BASIC;
      else if( _isRealLPLoaded )
      {
         assert(_simplifier == 0);
         assert(row < _solver.nRows());
         return _solver.getBasisRowStatus(row);
      }
      else
      {
         assert(row < _basisStatusRowsReal.size());
         return _basisStatusRowsReal[row];
      }
   }



   /// returns basis status for a single row
   SPxSolver::VarStatus SoPlex2::basisRowStatusReal(const SPxRowId& id) const
   {
      return basisRowStatusReal(idxReal(id));
   }



   /// returns basis status for a single column
   SPxSolver::VarStatus SoPlex2::basisColStatusReal(int col) const
   {
      assert(col >= 0);
      assert(col < numColsReal());

      // if no basis is available, return slack basis
      if( !_hasBasisReal )
      {
         if( lowerReal(col) > -realParam(SoPlex2::INFTY) )
            return SPxSolver::ON_LOWER;
         else if( upperReal(col) < realParam(SoPlex2::INFTY) )
            return SPxSolver::ON_UPPER;
         else
            return SPxSolver::ZERO;
      }
      else if( _isRealLPLoaded )
      {
         assert(_simplifier == 0);
         assert(col < _solver.nCols());
         return _solver.getBasisColStatus(col);
      }
      else
      {
         assert(col < _basisStatusColsReal.size());
         return _basisStatusColsReal[col];
      }
   }



   /// returns basis status for a single column
   SPxSolver::VarStatus SoPlex2::basisColStatusReal(const SPxColId& id) const
   {
      return basisColStatusReal(idxReal(id));
   }



   /// gets current basis
   void SoPlex2::getBasisReal(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[]) const
   {
      // if no basis is available, return slack basis
      if( !_hasBasisReal )
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
      else if( _isRealLPLoaded )
      {
         assert(_simplifier == 0);
         assert(numRowsReal() == _solver.nRows());
         assert(numColsReal() == _solver.nCols());

         (void)_solver.getBasis(rows, cols);
      }
      else
      {
         assert(numRowsReal() == _basisStatusRowsReal.size());
         assert(numColsReal() == _basisStatusColsReal.size());

         for( int i = numRowsReal() - 1; i >= 0; i-- )
            rows[i] = _basisStatusRowsReal[i];

         for( int i = numColsReal() - 1; i >= 0; i-- )
            cols[i] = _basisStatusColsReal[i];
      }
   }



   /// sets starting basis via arrays of statuses
   void SoPlex2::setBasisReal(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[])
   {
      if( _isRealLPLoaded )
      {
         assert(numRowsReal() == _solver.nRows());
         assert(numColsReal() == _solver.nCols());

         ///@todo check whether this has been successful and adjust _hasBasisReal accordingly
         _solver.setBasis(rows, cols);
      }
      else
      {
         _basisStatusRowsReal.reSize(numRowsReal());
         _basisStatusColsReal.reSize(numColsReal());

         for( int i = numRowsReal() - 1; i >= 0; i-- )
            _basisStatusRowsReal[i] = rows[i];

         for( int i = numColsReal() - 1; i >= 0; i-- )
            _basisStatusColsReal[i] = cols[i];

      }

      _hasBasisReal = true;
   }



   /// clears starting basis
   void SoPlex2::clearBasisReal()
   {
      if( _isRealLPLoaded )
         _solver.reLoad();

      _hasBasisReal = false;
   }



   /// is an advanced starting basis available?
   bool SoPlex2::hasBasisRational() const
   {
      return _hasBasisRational;
   }



   /// returns basis status for a single row
   SPxSolver::VarStatus SoPlex2::basisRowStatusRational(int row) const
   {
      assert(row >= 0);
      assert(row < numRowsRational());

      // if no basis is available, return slack basis
      if( !_hasBasisRational )
         return SPxSolver::BASIC;
      else
      {
         assert(row < _basisStatusRowsRational.size());
         return _basisStatusRowsRational[row];
      }
   }



   /// returns basis status for a single row
   SPxSolver::VarStatus SoPlex2::basisRowStatusRational(const SPxRowId& id) const
   {
      return basisRowStatusRational(idxRational(id));
   }



   /// returns basis status for a single column
   SPxSolver::VarStatus SoPlex2::basisColStatusRational(int col) const
   {
      assert(col >= 0);
      assert(col < numColsRational());

      // if no basis is available, return slack basis
      if( !_hasBasisRational )
      {
         if( lowerRational(col) > -realParam(SoPlex2::INFTY) )
            return SPxSolver::ON_LOWER;
         else if( upperRational(col) < realParam(SoPlex2::INFTY) )
            return SPxSolver::ON_UPPER;
         else
            return SPxSolver::ZERO;
      }
      else
      {
         assert(col < _basisStatusColsRational.size());
         return _basisStatusColsRational[col];
      }
   }



   /// returns basis status for a single column
   SPxSolver::VarStatus SoPlex2::basisColStatusRational(const SPxColId& id) const
   {
      return basisColStatusRational(idxRational(id));
   }



   /// gets current basis
   void SoPlex2::getBasisRational(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[]) const
   {
      // if no basis is available, return slack basis
      if( !_hasBasisRational )
      {
         for( int i = numRowsRational() - 1; i >= 0; i-- )
            rows[i] = SPxSolver::BASIC;

         for( int i = numColsRational() - 1; i >= 0; i-- )
         {
            if( lowerRational(i) > -realParam(SoPlex2::INFTY) )
               cols[i] = SPxSolver::ON_LOWER;
            else if( upperRational(i) < realParam(SoPlex2::INFTY) )
               cols[i] = SPxSolver::ON_UPPER;
            else
               cols[i] = SPxSolver::ZERO;
         }
      }
      else
      {
         assert(numRowsRational() == _basisStatusRowsRational.size());
         assert(numColsRational() == _basisStatusColsRational.size());

         for( int i = numRowsRational() - 1; i >= 0; i-- )
            rows[i] = _basisStatusRowsRational[i];

         for( int i = numColsRational() - 1; i >= 0; i-- )
            cols[i] = _basisStatusColsRational[i];
      }
   }



   /// sets starting basis via arrays of statuses
   void SoPlex2::setBasisRational(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[])
   {
      _basisStatusRowsRational.reSize(numRowsRational());
      _basisStatusColsRational.reSize(numColsRational());

      for( int i = numRowsRational() - 1; i >= 0; i-- )
         _basisStatusRowsRational[i] = rows[i];

      for( int i = numColsRational() - 1; i >= 0; i-- )
         _basisStatusColsRational[i] = cols[i];

      _hasBasisRational = true;
   }



   /// clears starting basis
   void SoPlex2::clearBasisRational()
   {
      _hasBasisRational = false;
   }


#if 0
   /// time spent in factorizations
   Real SoPlex2::factorTime() const
   {
   }



   /// number of factorizations performed
   int SoPlex2::factorCount() const
   {
   }



   /// time spent in solves
   Real SoPlex2::luSolveTime() const
   {
   }



   /// number of solves performed
   int SoPlex2::luSolveCount() const
   {
   }
#endif


   /// number of iterations since last call to solve
   int SoPlex2::numIterations() const
   {
      return _statistics->iterations;
   }



   ///
   Real SoPlex2::solveTime() const
   {
       return _solver.time();
   }



   /// statistical information in form of a string
   std::string SoPlex2::statisticString() const
   {
      return _solver.statistics();
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



   /// name of scaling method before simplifier
   const char* SoPlex2::getFirstScalerName()
   {
      if( _firstScaler )
         return _firstScaler->getName();
      else
         return "none";
   }



   /// name of scaling method after simplifier
   const char* SoPlex2::getSecondScalerName()
   {
      if( _secondScaler )
         return _secondScaler->getName();
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



   /// reads real LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool SoPlex2::readFileReal(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
   {
      assert(_realLP != 0);

      // clear statistics
      _statistics->clearAllData();

      // start timing
      _statistics->readingTime.start();

      // read
      bool success = _realLP->readFile(filename, rowNames, colNames, intVars);
      setIntParam(SoPlex2::OBJSENSE, (_realLP->spxSense() == SPxLPReal::MAXIMIZE ? SoPlex2::OBJSENSE_MAXIMIZE : SoPlex2::OBJSENSE_MINIMIZE), true, true);
      _hasBasisReal = false;

      // stop timing
      _statistics->readingTime.stop();

      return success;
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



   /// reads basis information from \p filename and returns true on success; if \p rowNames and \p colNames are \c NULL,
   /// default names are assumed; returns true on success
   bool SoPlex2::readBasisFileReal(const char* filename, const NameSet* rowNames, const NameSet* colNames)
   {
      assert(_realLP != 0);

      // start timing
      _statistics->readingTime.start();

      // read
      if( !_isRealLPLoaded )
      {
         assert(_realLP != &_solver);

         _solver.loadLP(*_realLP);
         spx_free(_realLP);
         _realLP = &_solver;
         _isRealLPLoaded = true;
      }
      _hasBasisReal = _solver.readBasisFile(filename, rowNames, colNames);

      // stop timing
      _statistics->readingTime.stop();

      return _hasBasisReal;
   }



   /// writes basis information to \p filename; if \p rowNames and \p colNames are \c NULL, default names are used;
   /// returns true on success
   bool SoPlex2::writeBasisFileReal(const char* filename, const NameSet* rowNames, const NameSet* colNames)
   {
      if( !_isRealLPLoaded )
      {
         assert(_realLP != &_solver);

         _solver.loadLP(*_realLP);
         spx_free(_realLP);
         _realLP = &_solver;
         _isRealLPLoaded = true;

         if( _hasBasisReal )
         {
            assert(_basisStatusRowsReal.size() == numRowsReal());
            assert(_basisStatusColsReal.size() == numColsReal());

            ///@todo this should not fail even if the basis is invalid (wrong dimension or wrong number of basic
            ///      entries); fix either in SPxSolver or in SPxBasis
            _solver.setBasis(_basisStatusRowsReal.get_const_ptr(), _basisStatusColsReal.get_const_ptr());
         }
      }

      return _solver.writeBasisFile(filename, rowNames, colNames);
   }



   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void SoPlex2::writeStateReal(const char* filename, const NameSet* rowNames, const NameSet* colNames)
   {
      writeFileReal(filename);
      writeBasisFileReal(filename);
      // @todo write settings file
   }



   /// reads rational LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool SoPlex2::readFileRational(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
   {
      assert(_rationalLP != 0);

      // clear statistics
      _statistics->clearAllData();

      // start timing
      _statistics->readingTime.start();

      // read
      bool success = _rationalLP->readFile(filename, rowNames, colNames, intVars);
      setIntParam(SoPlex2::OBJSENSE, (_rationalLP->spxSense() == SPxLPRational::MAXIMIZE ? SoPlex2::OBJSENSE_MAXIMIZE : SoPlex2::OBJSENSE_MINIMIZE), true, true);
      _hasBasisRational = false;

      // stop timing
      _statistics->readingTime.stop();

      return success;
   }



   /// writes rational LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer; returns true on success
   bool SoPlex2::writeFileRational(const char* filename, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars) const
   {
      ///@todo implement return value
      _rationalLP->writeFile(filename, rowNames, colNames, intVars);
      return true;
   }



   /// reads basis information from \p filename and returns true on success; if \p rowNames and \p colNames are \c NULL,
   /// default names are assumed; returns true on success
   bool SoPlex2::readBasisFileRational(const char* filename, const NameSet* rowNames, const NameSet* colNames)
   {
      assert(filename != 0);

      // start timing
      _statistics->readingTime.start();

      // read
      spxifstream file(filename);

      if( !file )
         return false;

      // get problem size
      int numRows = numRowsRational();
      int numCols = numColsRational();

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
            DataKey key = colIdRational(j);
            tmpColNames->add(key, name.str().c_str());
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
            DataKey key = rowIdRational(i);
            tmpRowNames->add(key, name.str().c_str());
         }

         rowNamesPtr = tmpRowNames;
      }

      // initialize with default slack basis
      _basisStatusRowsRational.reSize(numRows);
      _basisStatusColsRational.reSize(numCols);

      for( int i = 0; i < numRows; i++ )
         _basisStatusRowsRational[i] = SPxSolver::BASIC;

      for( int i = 0; i < numCols; i++ )
      {
         if( lowerRational(i) == upperRational(i) )
            _basisStatusColsRational[i] = SPxSolver::FIXED;
         else if( lowerRational(i) <= -realParam(SoPlex2::INFTY) && upperRational(i) >= realParam(SoPlex2::INFTY) )
            _basisStatusColsRational[i] = SPxSolver::ZERO;
         else if( lowerRational(i) <= -realParam(SoPlex2::INFTY) )
            _basisStatusColsRational[i] = SPxSolver::ON_UPPER;
         else
            _basisStatusColsRational[i] = SPxSolver::ON_LOWER;
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
               _basisStatusColsRational[c] = SPxSolver::BASIC;
               _basisStatusRowsRational[r] = (lhsRational(r) == rhsRational(r))
                  ? SPxSolver::FIXED
                  : SPxSolver::ON_UPPER;
            }
            else if( !strcmp(mps.field1(), "XL") )
            {
               _basisStatusColsRational[c] = SPxSolver::BASIC;
               _basisStatusRowsRational[r] = (lhsRational(r) == rhsRational(r))
                  ? SPxSolver::FIXED
                  : SPxSolver::ON_LOWER;
            }
            else if( !strcmp(mps.field1(), "UL") )
            {
               _basisStatusColsRational[c] = SPxSolver::ON_UPPER;
            }
            else if( !strcmp(mps.field1(), "LL") )
            {
               _basisStatusColsRational[c] = SPxSolver::ON_LOWER;
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

      _hasBasisRational = !mps.hasError();

      // stop timing
      _statistics->readingTime.stop();

      return _hasBasisRational;
   }



   /// writes basis information to \p filename; if \p rowNames and \p colNames are \c NULL, default names are used;
   /// returns true on success
   bool SoPlex2::writeBasisFileRational(const char* filename, const NameSet* rowNames, const NameSet* colNames)
   {
      assert(filename != 0);

      std::ofstream file(filename);
      if( file == 0 )
         return false;

      file.setf(std::ios::left);
      file << "NAME  " << filename << "\n";

      // do not write basis if there is none
      if( !_hasBasisRational )
      {
         file << "ENDATA\n";
         return true;
      }

      // start writing
      int numRows = _basisStatusRowsRational.size();
      int numCols = _basisStatusColsRational.size();
      int row = 0;

      for( int col = 0; col < numCols; col++ )
      {
         assert(_basisStatusColsRational[col] != SPxSolver::UNDEFINED);

         if( _basisStatusColsRational[col] == SPxSolver::BASIC )
         {
            // find nonbasic row
            for( ; row < numRows; row++ )
            {
               assert(_basisStatusRowsRational[row] != SPxSolver::UNDEFINED);
               if( _basisStatusRowsRational[row] != SPxSolver::BASIC )
                  break;
            }

            assert(row != numRows);

            file << (_basisStatusRowsRational[row] == SPxSolver::ON_UPPER ? " XU " : " XL ");

            file << std::setw(8);
            if( colNames != 0 && colNames->has(colIdRational(col)) )
               file << (*colNames)[colIdRational(col)];
            else
               file << "x" << col;

            file << "       ";
            if( rowNames != 0 && rowNames->has(rowIdRational(row)) )
               file << (*rowNames)[rowIdRational(row)];
            else
               file << "C" << row;

            file << "\n";
            row++;
         }
         else
         {
            if( _basisStatusColsRational[col] == SPxSolver::ON_UPPER )
            {
               file << " UL ";

               file << std::setw(8);
               if( colNames != 0 && colNames->has(colIdRational(col)) )
                  file << (*colNames)[colIdRational(col)];
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
         assert(_basisStatusRowsRational[row] == SPxSolver::BASIC);
      }
#endif

      return true;
   }



#if 0
   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void SoPlex2::writeStateRational(const char* filename, const NameSet* rowNames, const NameSet* colNames) const
   {
   }
#endif



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
         _rationalLP->changeSense(value == SoPlex2::OBJSENSE_MAXIMIZE ? SPxLPRational::MAXIMIZE : SPxLPRational::MINIMIZE);
         _invalidateSolutionReal();
         _invalidateSolutionRational();
         break;

      // type of computational form, i.e., column or row representation
      case SoPlex2::REPRESENTATION:
         if( value != SoPlex2::REPRESENTATION_COLUMN && value != SoPlex2::REPRESENTATION_ROW )
            return false;
         _solver.setRep(value == SoPlex2::REPRESENTATION_COLUMN ? SPxSolver::COLUMN : SPxSolver::ROW);
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
         {
            _solver.setTerminationIter(value);
            break;
         }

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

      // type of scaler applied before simplification
      case SoPlex2::SCALER_BEFORE_SIMPLIFIER:
         switch( value )
         {
         case SCALER_OFF:
            _firstScaler = 0;
            break;
         case SCALER_UNIEQUI:
            _firstScaler = &_scalerUniequi;
            break;
         case SCALER_BIEQUI:
            _firstScaler = &_scalerBiequi;
            break;
         case SCALER_GEO1:
            _firstScaler = &_scalerGeo1;
            break;
         case SCALER_GEO8:
            _firstScaler = &_scalerGeo8;
            break;
         default:
            return false;
         }
         break;

      // type of scaler applied after simplification
      case SoPlex2::SCALER_AFTER_SIMPLIFIER:
         switch( value )
         {
         case SCALER_OFF:
            _secondScaler = 0;
            break;
         case SCALER_UNIEQUI:
            _secondScaler = &_scalerUniequi;
            break;
         case SCALER_BIEQUI:
            _secondScaler = &_scalerBiequi;
            break;
         case SCALER_GEO1:
            _secondScaler = &_scalerGeo1;
            break;
         case SCALER_GEO8:
            _secondScaler = &_scalerGeo8;
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
         _solver.setTerminationIter(value);
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
         _solver.setFeastol(value);
         break;

      // dual feasibility tolerance
      case SoPlex2::OPTTOL:
         _solver.setOpttol(value);
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
   bool SoPlex2::saveSettingsFile(const char* filename) const
   {
      assert(filename != 0);

      std::ofstream file(filename);
      if( file == 0 )
         return false;

      file.setf(std::ios::left);
      file << "# SoPlex version " << SOPLEX_VERSION / 100 << "." << (SOPLEX_VERSION / 10) % 10 << "." << SOPLEX_VERSION % 10 << "." << SOPLEX_SUBVERSION << "\n";

      for( int i = 0; i < SoPlex2::BOOLPARAM_COUNT; i++ )
      {
         file << "\n";
         file << "# " << _currentSettings->_boolParamDescription[i] << "\n";
         file << "# range {true, false}, default " << _currentSettings->_boolParamDefault[i] << "\n";
         file << "bool:" << _currentSettings->_boolParamName[i] << " = " << (_currentSettings->_boolParamValues[i] ? "true\n" : "false\n");
      }

      for( int i = 0; i < SoPlex2::INTPARAM_COUNT; i++ )
      {
         file << "\n";
         file << "# " << _currentSettings->_intParamDescription[i] << "\n";
         file << "# range [-2147483648,2147483647], default " << _currentSettings->_intParamDefault[i] << "\n";
         file << "int:" << _currentSettings->_intParamName[i] << " = " << _currentSettings->_intParamValues[i] << "\n";
      }

      for( int i = 0; i < SoPlex2::REALPARAM_COUNT; i++ )
      {
         file << "\n";
         file << "# " << _currentSettings->_realParamDescription[i] << "\n";
         file << "# range [" << _currentSettings->_realParamLower[i] << "," << _currentSettings->_realParamLower[i]
            << "], default " << _currentSettings->_realParamDefault[i] << "\n";
         file << "real:" << _currentSettings->_realParamName[i] << " = " << _currentSettings->_realParamValues[i] << "\n";
      }

      for( int i = 0; i < SoPlex2::RATIONALPARAM_COUNT; i++ )
      {
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



   /// prints statistics on real solution
   void SoPlex2::printSolutionStatisticsReal(std::ostream& os)
   {
      os << "Solution           : \n"
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

      printStatus(os, _statusReal);

      os << "Rational LP        : \n"
         << "  Objective sense  : " << (intParam(SoPlex2::OBJSENSE) == SoPlex2::OBJSENSE_MINIMIZE ? "minimize\n" : "maximize\n");
      _realLP->printProblemStatistics(os);

      printSolutionStatisticsReal(os);

      printSolvingStatistics(os);
   }



   /// prints complete rational statistics
   void SoPlex2::printStatisticsRational(std::ostream& os)
   {
      os << std::setprecision(2);

      printStatus(os, _statusRational);

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



   /// creates a permutation for removing rows/columns from an array of IDs
   void SoPlex2::_idToPerm(SPxId* id, int idSize, int* perm, int permSize) const
   {
      assert(id != 0);
      assert(idSize >= 0);
      assert(perm != 0);
      assert(permSize >= 0);

      for( int i = 0; i < permSize; i++ )
         perm[i] = i;

      for( int i = 0; i < idSize; i++ )
      {
         assert(idxReal(id[i]) >= 0);
         assert(idxReal(id[i]) < permSize);
         perm[idxReal(id[i])] = -1;
      }
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
      assert(_rationalLP != 0);

      assert(_realLP != &_solver || _isRealLPLoaded);
      assert(_realLP == &_solver || !_isRealLPLoaded);

      assert(!_hasBasisReal || _isRealLPLoaded || _basisStatusRowsReal.size() == numRowsReal());
      assert(!_hasBasisReal || _isRealLPLoaded || _basisStatusColsReal.size() == numColsReal());
      assert(!_hasBasisRational || _basisStatusRowsRational.size() == numRowsRational());
      assert(!_hasBasisRational || _basisStatusColsRational.size() == numColsRational());

#if 0
      // this is not required since within _solveRational() we currently transform to minimization
      assert(intParam(SoPlex2::OBJSENSE) != SoPlex2::OBJSENSE_MAXIMIZE || _realLP->spxSense() == SPxLPReal::MAXIMIZE);
      assert(intParam(SoPlex2::OBJSENSE) != SoPlex2::OBJSENSE_MINIMIZE || _realLP->spxSense() == SPxLPReal::MINIMIZE);
      assert(intParam(SoPlex2::OBJSENSE) != SoPlex2::OBJSENSE_MAXIMIZE || _rationalLP->spxSense() == SPxLPRational::MAXIMIZE);
      assert(intParam(SoPlex2::OBJSENSE) != SoPlex2::OBJSENSE_MINIMIZE || _rationalLP->spxSense() == SPxLPRational::MINIMIZE);
#endif

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



   /// parses one line in a settings file; returns true on success
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

               if( sscanf(paramValueString, "%"REAL_FORMAT, &value) == 1 && setRealParam((SoPlex2::RealParam)param, value) )
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



   /// invalidates real solution
   void SoPlex2::_invalidateSolutionReal()
   {
      _statusReal = SPxSolver::UNKNOWN;
      _hasPrimalReal = false;
      _hasPrimalrayReal = false;
      _hasDualReal = false;
      _hasDualfarkasReal = false;
   }



   /// invalidates rational solution
   void SoPlex2::_invalidateSolutionRational()
   {
      _statusRational = SPxSolver::UNKNOWN;
      _solRational._invalidate();
   }



   /// enables simplifier and scalers according to current parameters
   void SoPlex2::_enableSimplifierAndScalers()
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

      // type of scaler applied before simplification
      switch( intParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER) )
      {
      case SCALER_OFF:
         _firstScaler = 0;
         break;
      case SCALER_UNIEQUI:
         _firstScaler = &_scalerUniequi;
         break;
      case SCALER_BIEQUI:
         _firstScaler = &_scalerBiequi;
         break;
      case SCALER_GEO1:
         _firstScaler = &_scalerGeo1;
         break;
      case SCALER_GEO8:
         _firstScaler = &_scalerGeo8;
         break;
      default:
         break;
      }

      // type of scaler applied after simplification
      switch( intParam(SoPlex2::SCALER_AFTER_SIMPLIFIER) )
      {
      case SCALER_OFF:
         _secondScaler = 0;
         break;
      case SCALER_UNIEQUI:
         _secondScaler = &_scalerUniequi;
         break;
      case SCALER_BIEQUI:
         _secondScaler = &_scalerBiequi;
         break;
      case SCALER_GEO1:
         _secondScaler = &_scalerGeo1;
         break;
      case SCALER_GEO8:
         _secondScaler = &_scalerGeo8;
         break;
      default:
         break;
      }
   }



   /// disables simplifier and scalers
   void SoPlex2::_disableSimplifierAndScalers()
   {
      _simplifier = 0;
      _firstScaler = 0;
      _secondScaler = 0;
   }



   /// synchronizes rational solution with real solution
   void SoPlex2::_syncRationalSolution(bool snycPrimal, bool syncDual, bool syncBasis)
   {
      DVectorReal buffer;
      _solRational._invalidate();

      if( hasPrimalReal() )
      {
         _solRational._hasPrimal = true;

         buffer.reDim(numColsReal());
         getPrimalReal(buffer);
         _solRational._primal = buffer;

         buffer.reDim(numRowsReal());
         getSlacksReal(buffer);
         _solRational._slacks = buffer;
      }

      if( hasPrimalrayReal() )
      {
         _solRational._hasPrimalray = true;

         buffer.reDim(numColsReal());
         getPrimalrayReal(buffer);
         _solRational._primalray = buffer;
      }

      if( hasDualReal() )
      {
         _solRational._hasDual = true;

         buffer.reDim(numRowsReal());
         getDualReal(buffer);
         _solRational._dual = buffer;

         buffer.reDim(numColsReal());
         getRedcostReal(buffer);
         _solRational._redcost = buffer;
      }

      if( hasDualfarkasReal() )
      {
         _solRational._hasDualfarkas = true;

         buffer.reDim(numRowsReal());
         getDualfarkasReal(buffer);
         _solRational._dualfarkas = buffer;
      }

      if( hasBasisReal() )
      {
         _hasBasisRational = true;

         _basisStatusRowsRational.reSize(numRowsReal());
         _basisStatusColsRational.reSize(numColsReal());
         getBasisReal(_basisStatusRowsRational.get_ptr(), _basisStatusColsRational.get_ptr());
      }
      else
         _hasBasisRational = false;
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
