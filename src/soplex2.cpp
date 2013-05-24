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

namespace soplex
{
   /// class of parameter settings
   class SoPlex2::Settings
   {
   public:
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
            _boolParamDefault[SoPlex2::PARTIAL_PRICING] = false;

            // objective sense
            _intParamDefault[SoPlex2::OBJSENSE] = SoPlex2::OBJSENSE_MAXIMIZE;

            // type of computational form, i.e., column or row representation
            _intParamDefault[SoPlex2::REPRESENTATION] = SoPlex2::REPRESENTATION_COLUMN;

            // type of algorithm, i.e., enter or leave
            _intParamDefault[SoPlex2::ALGORITHM] = SoPlex2::ALGORITHM_LEAVE;

            // type of LU update
            _intParamDefault[SoPlex2::FACTOR_UPDATE_TYPE] = SoPlex2::FACTOR_UPDATE_TYPE_FT;

            ///@todo which value?
            // maximum number of updates before fresh factorization
            _intParamDefault[SoPlex2::FACTOR_UPDATE_MAX] = 200;

            // iteration limit (-1 if unlimited)
            _intParamDefault[SoPlex2::ITERLIMIT] = -1;

            // display frequency
            _intParamDefault[SoPlex2::DISPLAY_FREQ] = 100;

            // verbosity level
            _intParamDefault[SoPlex2::VERBOSITY] = SoPlex2::VERBOSITY_NORMAL;

            // type of simplifier
            _intParamDefault[SoPlex2::SIMPLIFIER] = SoPlex2::SIMPLIFIER_AUTO;

            // type of scaler applied before simplification
            _intParamDefault[SoPlex2::SCALER_BEFORE_SIMPLIFIER] = SoPlex2::SCALER_OFF;

            // type of scaler applied after simplification
            _intParamDefault[SoPlex2::SCALER_AFTER_SIMPLIFIER] = SoPlex2::SCALER_BIEQUI;

            // type of starter used to create crash basis
            _intParamDefault[SoPlex2::STARTER] = SoPlex2::STARTER_OFF;

            // type of pricer
            _intParamDefault[SoPlex2::PRICER] = SoPlex2::PRICER_QUICKSTEEP;

            // type of ratio test
            _intParamDefault[SoPlex2::RATIOTESTER] = SoPlex2::RATIOTESTER_FAST;

            ///@todo define suitable values depending on Real type
            // general zero tolerance
            _realParamLower[SoPlex2::EPSILON_ZERO] = DEFAULT_EPS_ZERO;
            _realParamUpper[SoPlex2::EPSILON_ZERO] = DEFAULT_EPS_ZERO;
            _realParamDefault[SoPlex2::EPSILON_ZERO] = DEFAULT_EPS_ZERO;

            ///@todo define suitable values depending on Real type
            // zero tolerance used in factorization
            _realParamLower[SoPlex2::EPSILON_FACTORIZATION] = DEFAULT_EPS_FACTOR;
            _realParamUpper[SoPlex2::EPSILON_FACTORIZATION] = DEFAULT_EPS_FACTOR;
            _realParamDefault[SoPlex2::EPSILON_FACTORIZATION] = DEFAULT_EPS_FACTOR;

            ///@todo define suitable values depending on Real type
            // zero tolerance used in factorization update
            _realParamLower[SoPlex2::EPSILON_UPDATE] = DEFAULT_EPS_UPDATE;
            _realParamUpper[SoPlex2::EPSILON_UPDATE] = DEFAULT_EPS_UPDATE;
            _realParamDefault[SoPlex2::EPSILON_UPDATE] = DEFAULT_EPS_UPDATE;

            ///@todo define suitable values depending on Real type
            // infinity threshold
            _realParamLower[SoPlex2::INFTY] = DEFAULT_INFINITY;
            _realParamUpper[SoPlex2::INFTY] = DEFAULT_INFINITY;
            _realParamDefault[SoPlex2::INFTY] = DEFAULT_INFINITY;

            // time limit in seconds (INFTY if unlimited)
            _realParamLower[SoPlex2::TIMELIMIT] = 0.0;
            _realParamUpper[SoPlex2::TIMELIMIT] = _realParamLower[SoPlex2::INFTY];
            _realParamDefault[SoPlex2::TIMELIMIT] = _realParamLower[SoPlex2::INFTY];

            // lower limit on objective value
            _realParamLower[SoPlex2::OBJLIMIT_LOWER] = -_realParamLower[SoPlex2::INFTY];
            _realParamUpper[SoPlex2::OBJLIMIT_LOWER] = _realParamLower[SoPlex2::INFTY];
            _realParamDefault[SoPlex2::OBJLIMIT_LOWER] = -_realParamLower[SoPlex2::INFTY];

            // upper limit on objective value
            _realParamLower[SoPlex2::OBJLIMIT_UPPER] = -_realParamLower[SoPlex2::INFTY];
            _realParamUpper[SoPlex2::OBJLIMIT_UPPER] = _realParamLower[SoPlex2::INFTY];
            _realParamDefault[SoPlex2::OBJLIMIT_UPPER] = _realParamLower[SoPlex2::INFTY];

            ///@todo define suitable values depending on Real type
            // threshold for activating iterative refinement
            _realParamLower[SoPlex2::IRTHRESHOLD] = 0.0;
            _realParamUpper[SoPlex2::IRTHRESHOLD] = 1.0;
            _realParamDefault[SoPlex2::IRTHRESHOLD] = 1e-12;

            // primal feasibility tolerance
            _rationalParamLower[SoPlex2::FEASTOL] = 0.0;
            _rationalParamUpper[SoPlex2::FEASTOL] = 1.0;
            _rationalParamDefault[SoPlex2::FEASTOL] = 1e-6;

            // dual feasibility tolerance
            _rationalParamLower[SoPlex2::OPTTOL] = 0.0;
            _rationalParamUpper[SoPlex2::OPTTOL] = 1.0;
            _rationalParamDefault[SoPlex2::OPTTOL] = 1e-6;

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
      }
   };

   bool SoPlex2::Settings::_defaultsAndBoundsInitialized = false;
   bool SoPlex2::Settings::_boolParamDefault[SoPlex2::BOOLPARAM_COUNT];
   int SoPlex2::Settings::_intParamDefault[SoPlex2::INTPARAM_COUNT];
   Real SoPlex2::Settings::_realParamLower[SoPlex2::REALPARAM_COUNT];
   Real SoPlex2::Settings::_realParamUpper[SoPlex2::REALPARAM_COUNT];
   Real SoPlex2::Settings::_realParamDefault[SoPlex2::REALPARAM_COUNT];
   Rational SoPlex2::Settings::_rationalParamLower[SoPlex2::RATIONALPARAM_COUNT];
   Rational SoPlex2::Settings::_rationalParamUpper[SoPlex2::RATIONALPARAM_COUNT];
   Rational SoPlex2::Settings::_rationalParamDefault[SoPlex2::RATIONALPARAM_COUNT];



   /// default constructor
   SoPlex2::SoPlex2()
      : _scalerUniequi(false)
      , _scalerBiequi(true)
      , _scalerGeo1(1)
      , _scalerGeo8(8)
      , _simplifier(0)
      , _firstScaler(0)
      , _secondScaler(0)
      , _starter(0)
      , _hasBasisReal(false)
   {
      // give lu factorization to solver
      _solver.setSolver(&_slufactor);

      // the real LP is initially stored in the solver
      _realLP = &_solver;
      _isRealLPLoaded = true;

      // construct the rational LP
      spx_alloc(_rationalLP);
      _rationalLP = new (_rationalLP) SPxLPRational();

      // initialize parameter settings to default
      spx_alloc(_currentSettings);
      _currentSettings = new (_currentSettings) Settings();
      setSettings(*_currentSettings, true);

      assert(_isConsistent());
   }



   /// assignment operator
   SoPlex2& SoPlex2::operator=(const SoPlex2& rhs)
   {
      assert(rhs._rationalLP != 0);

      if( this != &rhs )
      {
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

         // initialize pointers for simplifier, scalers, and starter
         setIntParam(SoPlex2::SIMPLIFIER, intParam(SoPlex2::SIMPLIFIER), true);
         setIntParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER, intParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER), true);
         setIntParam(SoPlex2::SCALER_AFTER_SIMPLIFIER, intParam(SoPlex2::SCALER_AFTER_SIMPLIFIER), true);
         setIntParam(SoPlex2::STARTER, intParam(SoPlex2::STARTER), true);

         // copy real LP if different from the LP in the solver
         if( rhs._realLP != &(rhs._solver) )
         {
            spx_alloc(_realLP);
            _realLP = new (_realLP) SPxLPReal(*(rhs._realLP));
         }
         else
            _realLP = &_solver;

         // copy rational LP
         spx_alloc(_rationalLP);
         _rationalLP = new (_rationalLP) SPxLPRational(*rhs._rationalLP);

         // copy boolean flags
         _isRealLPLoaded = rhs._isRealLPLoaded;
         _hasBasisReal = rhs._hasBasisReal;
      }
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
   typename LPRowReal::Type SoPlex2::rowTypeReal(int i) const
   {
      assert(_realLP != 0);
      return _realLP->rowType(i);
   }



   /// returns inequality type of row with identifier \p id
   typename LPRowReal::Type SoPlex2::rowTypeReal(const SPxRowId& id) const
   {
      assert(_realLP != 0);
      return _realLP->rowType(id);
   }



   /// gets column \p i
   void SoPlex2::getColReal(int i, LPCol& lpcol) const
   {
      assert(_realLP != 0);
      return _realLP->getCol(i, lpcol);
   }



   /// gets column with identifier \p id.
   void SoPlex2::getColReal(const SPxColId& id, LPCol& lpcol) const
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



   /// adds a single row
   void SoPlex2::addRowReal(const LPRowReal& lprow)
   {
      SPxRowId id;
      addRowReal(id, lprow);
   }



   /// adds a single row and gets its \p id
   void SoPlex2::addRowReal(SPxRowId& id, const LPRowReal& lprow)
   {
      assert(_realLP != 0);
      _realLP->addRow(id, lprow);

      if( _hasBasisReal && !_isRealLPLoaded )
         _basisStatusRowsReal.append(SPxSolver::BASIC);
   }



   /// adds multiple rows
   void SoPlex2::addRowsReal(const LPRowSetReal& lprowset)
   {
      assert(_realLP != 0);
      _realLP->addRows(lprowset);

      if( _hasBasisReal && !_isRealLPLoaded )
         _basisStatusRowsReal.append(lprowset.num(), SPxSolver::BASIC);
   }



   /// adds multiple rows and gets an array of their \p id 's
   void SoPlex2::addRowsReal(SPxRowId id[], const LPRowSetReal& lprowset)
   {
      assert(_realLP != 0);
      _realLP->addRows(id, lprowset);

      if( _hasBasisReal && !_isRealLPLoaded )
         _basisStatusRowsReal.append(lprowset.num(), SPxSolver::BASIC);
   }



   /// adds a single column
   void SoPlex2::addColReal(const LPCol& lpcol)
   {
      SPxColId id;
      addColReal(id, lpcol);
   }



   /// adds a single column and gets its \p id
   void SoPlex2::addColReal(SPxColId& id, const LPCol& lpcol)
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
   }



   /// replaces row with identifier \p id with \p lprow
   void SoPlex2::changeRowReal(SPxRowId id, const LPRowReal& lprow)
   {
      changeRowReal(idxReal(id), lprow);
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
   }



   /// changes left-hand side of row \p i to \p lhs
   void SoPlex2::changeLhsReal(int i, Real lhs)
   {
      assert(_realLP != 0);
      _realLP->changeLhs(i, lhs);

      if( _hasBasisReal && !_isRealLPLoaded && _basisStatusRowsReal[i] == SPxSolver::ON_LOWER && lhs <= -realParam(SoPlex2::INFTY) )
         _basisStatusRowsReal[i] = (rhsReal(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
   }



   /// changes left-hand side of row with identifier \p id to \p lhs
   void SoPlex2::changeLhsReal(SPxRowId id, Real lhs)
   {
      changeLhsReal(idxReal(id), lhs);
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
   }



   /// changes right-hand side of row \p i to \p rhs
   void SoPlex2::changeRhsReal(int i, Real rhs)
   {
      assert(_realLP != 0);
      _realLP->changeRhs(i, rhs);

      if( _hasBasisReal && !_isRealLPLoaded && _basisStatusRowsReal[i] == SPxSolver::ON_UPPER && rhs >= realParam(SoPlex2::INFTY) )
         _basisStatusRowsReal[i] = (lhsReal(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
   }



   /// changes right-hand of row with identifier \p id to \p rhs
   void SoPlex2::changeRhsReal(SPxRowId id, Real rhs)
   {
      changeRhsReal(idxReal(id), rhs);
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
   }



   /// changes left- and right-hand side of row with identifier \p id
   void SoPlex2::changeRangeReal(SPxRowId id, Real lhs, Real rhs)
   {
      changeRangeReal(idxReal(id), lhs, rhs);
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
   }



   /// replaces column with identifier \p id with \p lpcol
   void SoPlex2::changeColReal(SPxColId id, const LPColReal& lpcol)
   {
      changeColReal(idxReal(id), lpcol);
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
   }



   /// changes lower bound of column i to \p lower
   void SoPlex2::changeLowerReal(int i, Real lower)
   {
      assert(_realLP != 0);
      _realLP->changeLower(i, lower);

      if( _hasBasisReal && !_isRealLPLoaded && _basisStatusColsReal[i] == SPxSolver::ON_LOWER && lower <= -realParam(SoPlex2::INFTY) )
         _basisStatusColsReal[i] = (upperReal(i) < realParam(SoPlex2::INFTY)) ? SPxSolver::ON_UPPER : SPxSolver::ZERO;
   }



   /// changes lower bound of column with identifier \p id to \p lower
   void SoPlex2::changeLowerReal(SPxColId id, Real lower)
   {
      changeLowerReal(idxReal(id), lower);
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
   }



   /// changes \p i 'th upper bound to \p upper
   void SoPlex2::changeUpperReal(int i, Real upper)
   {
      assert(_realLP != 0);
      _realLP->changeUpper(i, upper);

      if( _hasBasisReal && !_isRealLPLoaded &&  _basisStatusColsReal[i] == SPxSolver::ON_UPPER && upper >= realParam(SoPlex2::INFTY) )
         _basisStatusColsReal[i] = (lowerReal(i) > -realParam(SoPlex2::INFTY)) ? SPxSolver::ON_LOWER : SPxSolver::ZERO;
   }



   /// changes upper bound of column with identifier \p id to \p upper
   void SoPlex2::changeUpperReal(SPxColId id, Real upper)
   {
      changeUpperReal(idxReal(id), upper);
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
   }



   /// changes bounds of column with identifier \p id to \p lower and \p upper
   void SoPlex2::changeBoundsReal(SPxColId id, Real lower, Real upper)
   {
      changeBoundsReal(idxReal(id), lower, upper);
   }



   /// changes objective function vector to \p obj
   void SoPlex2::changeObjReal(const VectorReal& obj)
   {
      assert(_realLP != 0);
      _realLP->changeObj(obj);
   }



   /// changes objective coefficient of column i to \p obj
   void SoPlex2::changeObjReal(int i, Real obj)
   {
      assert(_realLP != 0);
      _realLP->changeObj(i, obj);
   }



   /// changes objective coefficient of column with identifier \p id to \p obj
   void SoPlex2::changeObjReal(SPxColId id, Real obj)
   {
      assert(_realLP != 0);
      _realLP->changeObj(id, obj);
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
   }



   /// changes matrix entry identified by (\p rowid, \p colid) to \p val
   void SoPlex2::changeElementReal(SPxRowId rowid, SPxColId colid, Real val)
   {
      changeElementReal(idxReal(rowid), idxReal(colid), val);
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
   }



   /// removes row with identifier \p id
   void SoPlex2::removeRowReal(SPxRowId id)
   {
      removeRowReal(idxReal(id));
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
   }



   /// remove all rows with identifier in array \p id of size \p n; an array \p perm of size #numRowsReal() may be
   /// passed as buffer memory
   void SoPlex2::removeRowsReal(SPxRowId id[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numRowsReal());
         _idToPerm((SPxId*)id, n, p.get_ptr(), numRowsReal());
         removeRowsReal(p.get_ptr());
      }
      else
      {
         _idToPerm((SPxId*)id, n, perm, numRowsReal());
         removeRowsReal(perm);
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
         removeRowsReal(p.get_ptr());
      }
      else
      {
         _idxToPerm(idx, n, perm, numRowsReal());
         removeRowsReal(perm);
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
         removeRowsReal(p.get_ptr());
      }
      else
      {
         _rangeToPerm(start, end, perm, numRowsReal());
         removeRowsReal(perm);
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
   }



   /// removes column with identifier \p id
   void SoPlex2::removeColReal(SPxColId id)
   {
      removeColReal(idxReal(id));
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
   }



   /// remove all columns with identifier in array \p id of size \p n; an array \p perm of size #numColsReal() may be
   /// passed as buffer memory
   void SoPlex2::removeColsReal(SPxColId id[], int n, int perm[])
   {
      if( perm == 0 )
      {
         DataArray< int > p(numColsReal());
         _idToPerm((SPxId*)id, n, p.get_ptr(), numColsReal());
         removeColsReal(p.get_ptr());
      }
      else
      {
         _idToPerm((SPxId*)id, n, perm, numColsReal());
         removeColsReal(perm);
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
         removeColsReal(p.get_ptr());
      }
      else
      {
         _idxToPerm(idx, n, perm, numColsReal());
         removeColsReal(perm);
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
         removeColsReal(p.get_ptr());
      }
      else
      {
         _rangeToPerm(start, end, perm, numColsReal());
         removeColsReal(perm);
      }
   }



   /// clears the LP
   void SoPlex2::clearLPReal()
   {
      assert(_realLP != 0);

      _realLP->clear();
      _hasBasisReal = false;
   }



   /// solves real LP
   SPxSolver::Status SoPlex2::solveReal()
   {
      assert(_isConsistent());

      // will preprocessing be applied? (only if no basis is available)
      bool applyPreprocessing = (_firstScaler != 0 || _simplifier != 0 || _secondScaler != 0) && !_hasBasisReal;

      if( _isRealLPLoaded )
      {
         assert(_realLP == &_solver);
         assert(applyPreprocessing || _hasBasisReal);

         if( applyPreprocessing )
         {
            // why is this done? _realLP is now a copy of the loaded LP of _solver
            spx_alloc(_realLP);
            _realLP = new (_realLP) SPxLPReal(_solver);
            _isRealLPLoaded = false; // why do we set this to false?
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
         assert(applyPreprocessing || _hasBasisReal);

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
         ///@todo this should become a separate method that implements starter and auto pricing
         _solver.solve();
      }

      // evaluate result
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

      return statusReal();
   }



   /// returns the current status
   SPxSolver::Status SoPlex2::statusReal() const
   {
      if( _simplifier == 0 || _simplifier->result() == SPxSimplifier::OKAY )
         return _solver.status();
      else if( _simplifier->result() == SPxSimplifier::INFEASIBLE )
         return SPxSolver::INFEASIBLE;
      else if( _simplifier->result() == SPxSimplifier::DUAL_INFEASIBLE )
         return SPxSolver::INForUNBD;
      else if( _simplifier->result() == SPxSimplifier::UNBOUNDED )
         return SPxSolver::UNBOUNDED;
      else if( _simplifier->result() == SPxSimplifier::VANISHED )
         return SPxSolver::OPTIMAL;
      else
         return SPxSolver::UNKNOWN;
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



   /// is a primal feasible solution available?
   bool SoPlex2::hasPrimalReal() const
   {
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

      return statusReal() == SPxSolver::OPTIMAL || (_simplifier == 0 &&
         (_solver.basis().status() == SPxBasis::PRIMAL || _solver.basis().status() == SPxBasis::UNBOUNDED) &&
         _solver.shift() < 10.0 * realParam(SoPlex2::EPSILON_ZERO));
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



   /// gets the primal ray if LP is unbounded; returns true on success
   bool SoPlex2::getPrimalrayReal(VectorReal& vector) const
   {
      if( statusReal() == SPxSolver::UNBOUNDED && _simplifier == 0 )
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

      return statusReal() == SPxSolver::OPTIMAL || (_simplifier == 0 &&
         (_solver.basis().status() == SPxBasis::DUAL || _solver.basis().status() == SPxBasis::INFEASIBLE) &&
         _solver.shift() < 10.0 * realParam(SoPlex2::EPSILON_ZERO));
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
            _solver.getDual(vector);

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



   /// gets the Farkas proof if LP is infeasible; returns true on success
   bool SoPlex2::getDualfarkasReal(VectorReal& vector) const
   {
      if( statusReal() == SPxSolver::INFEASIBLE && _simplifier == 0 )
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

      for( int i = numColsReal(); i >= 0; i-- )
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

      DVectorReal activity = _realLP->computePrimalActivity(primal);

      maxviol = 0.0;
      sumviol = 0.0;

      for( int i = numRowsReal(); i >= 0; i-- )
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



   /// gets current basis and returns solver status
   void SoPlex2::getBasisReal(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[]) const
   {
      // if no basis is available, return slack basis
      if( !_hasBasisReal )
      {
         for( int i = numRowsReal(); i >= 0; i-- )
            rows[i] = SPxSolver::BASIC;

         for( int i = numColsReal(); i >= 0; i-- )
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

         for( int i = numRowsReal(); i >= 0; i-- )
            rows[i] = _basisStatusRowsReal[i];

         for( int i = numColsReal(); i >= 0; i-- )
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

         for( int i = numRowsReal(); i >= 0; i-- )
            _basisStatusRowsReal[i] = rows[i];

         for( int i = numColsReal(); i >= 0; i-- )
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



   /// number of iterations since last call to solve
   int SoPlex2::numIterations() const
   {
   }



   /// statistical information in form of a string
   std::string SoPlex2::statisticString() const
   {
   }
#endif


   /// reads real LP in LP or MPS format from file and returns true on success; gets row names, column names, and
   /// integer variables if desired
   bool SoPlex2::readFileReal(const char* filename, NameSet* rowNames, NameSet* colNames, DIdxSet* intVars)
   {
      assert(_realLP != 0);
      bool success = _realLP->readFile(filename, rowNames, colNames, intVars);
      _hasBasisReal = false;
   }



   /// writes real LP to file; LP or MPS format is chosen from the extension in \p filename; if \p rowNames and \p
   /// colNames are \c NULL, default names are used; if \p intVars is not \c NULL, the variables contained in it are
   /// marked as integer
   void SoPlex2::writeFileReal(const char* filename, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars) const
   {
      return _realLP->writeFile(filename, rowNames, colNames, intVars);
   }



   /// reads basis information from \p filename and returns true on success; if \p rowNames and \p colNames are \c NULL,
   /// default names are assumed
   bool SoPlex2::readBasisFileReal(const char* filename, const NameSet* rowNames, const NameSet* colNames)
   {
      assert(_realLP != 0);

      if( !_isRealLPLoaded )
      {
         assert(_realLP != &_solver);

         _solver.loadLP(*_realLP);
         spx_free(_realLP);
         _realLP = &_solver;
         _isRealLPLoaded = true;
      }

      _hasBasisReal = _solver.readBasisFile(filename, rowNames, colNames);
   }



   /// writes basis information to \p filename; if \p rowNames and \p colNames are \c NULL, default names are used
   void SoPlex2::writeBasisFileReal(const char* filename, const NameSet* rowNames, const NameSet* colNames)
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

      _solver.writeBasisFile(filename, rowNames, colNames);
   }



#if 0
   /// writes internal LP, basis information, and parameter settings; if \p rowNames and \p colNames are \c NULL,
   /// default names are used
   void SoPlex2::writeStateReal(const char* filename, const NameSet* rowNames, const NameSet* colNames) const
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
   bool SoPlex2::setBoolParam(const BoolParam param, const bool value, bool quiet)
   {
      assert(param >= 0);
      assert(param < SoPlex2::BOOLPARAM_COUNT);
      assert(_isConsistent());

      if( value == boolParam(param) )
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
      assert(_isConsistent());
      return true;
   }



   /// sets integer parameter value; returns true on success
   bool SoPlex2::setIntParam(const IntParam param, const int value, bool quiet)
   {
      assert(param >= 0);
      assert(param < INTPARAM_COUNT);
      assert(_isConsistent());

      if( value == intParam(param) )
         return true;

      switch( param )
      {
      // objective sense
      case SoPlex2::OBJSENSE:
         assert(value == SoPlex2::OBJSENSE_MAXIMIZE || value == SoPlex2::OBJSENSE_MINIMIZE);
         _realLP->changeSense(value == SoPlex2::OBJSENSE_MAXIMIZE ? SPxLPReal::MAXIMIZE : SPxLPReal::MINIMIZE);
         if( _rationalLP != 0 )
            _rationalLP->changeSense(value == SoPlex2::OBJSENSE_MAXIMIZE ? SPxLPRational::MAXIMIZE : SPxLPRational::MINIMIZE);
         break;

      // type of computational form, i.e., column or row representation
      case SoPlex2::REPRESENTATION:
         assert(value == SoPlex2::REPRESENTATION_COLUMN || value == SoPlex2::REPRESENTATION_ROW);
         _solver.setRep(value == SoPlex2::REPRESENTATION_COLUMN ? SPxSolver::COLUMN : SPxSolver::ROW);
         break;

      // type of algorithm, i.e., enter or leave
      case SoPlex2::ALGORITHM:
         assert(value == SoPlex2::ALGORITHM_ENTER || value == SoPlex2::ALGORITHM_LEAVE);
         _solver.setType(value == SoPlex2::ALGORITHM_ENTER ? SPxSolver::ENTER : SPxSolver::LEAVE);
         break;

      // type of LU update
      case SoPlex2::FACTOR_UPDATE_TYPE:
         assert(value == SoPlex2::FACTOR_UPDATE_TYPE_ETA || value == SoPlex2::FACTOR_UPDATE_TYPE_FT);
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
      case SoPlex2::DISPLAY_FREQ:
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
            break;
         default:
            return false;
         }

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

      default:
         return false;
      }

      _currentSettings->_intParamValues[param] = value;
      assert(_isConsistent());
      return true;
   }



   /// sets real parameter value; returns true on success
   bool SoPlex2::setRealParam(const RealParam param, const Real value, bool quiet)
   {
      assert(param >= 0);
      assert(param < REALPARAM_COUNT);
      assert(value >= Settings::_realParamLower[param]);
      assert(value <= _currentSettings->_realParamUpper[param]);
      assert(_isConsistent());

      if( value == realParam(param) )
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

      // zero tolerance used in factorization update
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

      // lower limit on objective value
      case SoPlex2::OBJLIMIT_LOWER:
         ///@todo implement for both objective senses in SPxSolver::terminate()
         return false;

      // upper limit on objective value
      case SoPlex2::OBJLIMIT_UPPER:
         ///@todo implement for both objective senses in SPxSolver::terminate()
         if( intParam(SoPlex2::OBJSENSE) == SoPlex2::OBJSENSE_MINIMIZE )
         {
            _solver.setTerminationValue(value);
            break;
         }
         else
            return false;

      // threshold for activating iterative refinement
      case SoPlex2::IRTHRESHOLD:
         _solver.setIrthreshold(value);
         break;

      default:
         return false;
      }

      _currentSettings->_realParamValues[param] = value;
      assert(_isConsistent());
      return true;
   }



   /// sets rational parameter value; returns true on success
   bool SoPlex2::setRationalParam(const RationalParam param, const Rational value, bool quiet)
   {
      assert(param >= 0);
      assert(param < RATIONALPARAM_COUNT);
      assert(value >= _currentSettings->_rationalParamLower[param]);
      assert(value <= _currentSettings->_rationalParamUpper[param]);
      assert(_isConsistent());

      if( value == rationalParam(param) )
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
      assert(_isConsistent());
      return true;
   }



   /// sets parameter settings; returns true on success
   bool SoPlex2::setSettings(const Settings& settings, bool quiet)
   {
      *_currentSettings = settings;

      for( int i = 0; i < SoPlex2::BOOLPARAM_COUNT; i++ )
         setBoolParam((BoolParam)i, _currentSettings->_boolParamValues[i], quiet);

      for( int i = 0; i < SoPlex2::INTPARAM_COUNT; i++ )
         setIntParam((IntParam)i, _currentSettings->_intParamValues[i], quiet);

      for( int i = 0; i < SoPlex2::REALPARAM_COUNT; i++ )
         setRealParam((RealParam)i, _currentSettings->_realParamValues[i], quiet);

      for( int i = 0; i < SoPlex2::RATIONALPARAM_COUNT; i++ )
         setRationalParam((RationalParam)i, _currentSettings->_rationalParamValues[i], quiet);
   }



   /// checks consistency
   bool SoPlex2::_isConsistent() const
   {
      assert(_currentSettings != 0);
      assert(_realLP != 0);

      assert(_realLP != &_solver || _isRealLPLoaded);
      assert(_realLP == &_solver || !_isRealLPLoaded);

      assert(!_hasBasisReal || _isRealLPLoaded || _basisStatusRowsReal.size() == numRowsReal());
      assert(!_hasBasisReal || _isRealLPLoaded || _basisStatusColsReal.size() == numColsReal());

      assert(intParam(SoPlex2::OBJSENSE) != SoPlex2::OBJSENSE_MAXIMIZE || _realLP->spxSense() == SPxLPReal::MAXIMIZE);
      assert(intParam(SoPlex2::OBJSENSE) != SoPlex2::OBJSENSE_MINIMIZE || _realLP->spxSense() == SPxLPReal::MINIMIZE);
      assert(intParam(SoPlex2::OBJSENSE) != SoPlex2::OBJSENSE_MAXIMIZE || _rationalLP != 0 || _rationalLP->spxSense() == SPxLPRational::MAXIMIZE);
      assert(intParam(SoPlex2::OBJSENSE) != SoPlex2::OBJSENSE_MINIMIZE || _rationalLP != 0 || _rationalLP->spxSense() == SPxLPRational::MINIMIZE);

      assert(intParam(SoPlex2::SIMPLIFIER) != SoPlex2::SIMPLIFIER_OFF || _simplifier == 0);
      assert(intParam(SoPlex2::SIMPLIFIER) == SoPlex2::SIMPLIFIER_OFF || _simplifier != 0);

      assert(intParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER) != SoPlex2::SCALER_OFF || _firstScaler == 0);
      assert(intParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER) == SoPlex2::SCALER_OFF || _firstScaler != 0);

      assert(intParam(SoPlex2::SCALER_AFTER_SIMPLIFIER) != SoPlex2::SCALER_OFF || _secondScaler == 0);
      assert(intParam(SoPlex2::SCALER_AFTER_SIMPLIFIER) == SoPlex2::SCALER_OFF || _secondScaler != 0);
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
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
