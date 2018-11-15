/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  soplex.hpp
 * @brief General templated functions for SoPlex
 */

/// default setting for LU refactorization interval
#define DEFAULT_REFACTOR_INTERVAL 200


/// returns number of columns
template <class R>
int SoPlexBase<R>::numCols() const
{
  assert(_realLP != nullptr);
  return _realLP->nCols();
}

/// returns number of rows
template <class R>
int SoPlexBase<R>::numRows() const
{
  assert(_realLP != nullptr);
  return _realLP->nRows();
}

/// returns number of nonzeros
template <class R>
int SoPlexBase<R>::numNonzeros() const
{
  assert(_realLP != 0);
  return _realLP->nNzos();
}

/// gets the primal solution vector if available; returns true on success
template <class R>
bool SoPlexBase<R>::getPrimal(VectorBase<R>& vector)
{
  if( hasPrimal() && vector.dim() >= numCols() )
    {
      _syncRealSolution();
      _solReal.getPrimalSol(vector);
      return true;
    }
  else
    return false;
}

/// gets the primal ray if available; returns true on success
template <class R>
bool SoPlexBase<R>::getPrimalRay(VectorBase<R>& vector)
{
  if( hasPrimalRay() && vector.dim() >= numCols() )
    {
      _syncRealSolution();
      _solReal.getPrimalRaySol(vector);
      return true;
    }
  else
    return false;
}


/// gets the dual solution vector if available; returns true on success
template <class R>
bool SoPlexBase<R>::getDual(VectorBase<R>& vector)
{
  if( hasDual() && vector.dim() >= numRows() )
    {
      _syncRealSolution();
      _solReal.getDualSol(vector);
      return true;
    }
  else
    return false;
}



/// gets the Farkas proof if available; returns true on success
template <class R>
bool SoPlexBase<R>::getDualFarkas(VectorBase<R>& vector)
{
  if( hasDualFarkas() && vector.dim() >= numRows() )
    {
      _syncRealSolution();
      _solReal.getDualFarkasSol(vector);
      return true;
    }
  else
    return false;
}


/// gets violation of bounds; returns true on success
template <class R>
bool SoPlexBase<R>::getBoundViolation(R& maxviol, R& sumviol)
{
  if( !isPrimalFeasible() )
    return false;

  _syncRealSolution();
  VectorReal& primal = _solReal._primal;
  assert(primal.dim() == numCols());

  maxviol = 0.0;
  sumviol = 0.0;

  for( int i = numCols() - 1; i >= 0; i-- )
    {
      Real lower = _realLP->lowerUnscaled(i);
      Real upper = _realLP->upperUnscaled(i);
      Real viol = lower - primal[i];
      if( viol > 0.0 )
        {
          sumviol += viol;
          if( viol > maxviol )
            maxviol = viol;
        }

      viol = primal[i] - upper;
      if( viol > 0.0 )
        {
          sumviol += viol;
          if( viol > maxviol )
            maxviol = viol;
        }
    }

  return true;
}

/// gets the vector of reduced cost values if available; returns true on success
template <class R>
bool SoPlexBase<R>::getRedCost(VectorBase<R>& vector)
{
  if( hasDual() && vector.dim() >= numCols() )
    {
      _syncRealSolution();
      _solReal.getRedCostSol(vector);
      return true;
    }
  else
    return false;
}

/// gets violation of constraints; returns true on success
template <class R>
bool SoPlexBase<R>::getRowViolation(R& maxviol, R& sumviol)
{
  if( !isPrimalFeasible() )
    return false;

  _syncRealSolution();
  VectorReal& primal = _solReal._primal;
  assert(primal.dim() == numCols());

  DVectorReal activity(numRows());
  _realLP->computePrimalActivity(primal, activity, true);
  maxviol = 0.0;
  sumviol = 0.0;

  for( int i = numRows() - 1; i >= 0; i-- )
    {
      Real lhs = _realLP->lhsUnscaled(i);
      Real rhs = _realLP->rhsUnscaled(i);

      Real viol = lhs - activity[i];
      if( viol > 0.0 )
        {
          sumviol += viol;
          if( viol > maxviol )
            maxviol = viol;
        }

      viol = activity[i] - rhs;
      if( viol > 0.0 )
        {
          sumviol += viol;
          if( viol > maxviol )
            maxviol = viol;
        }
    }

  return true;
}

/// gets violation of dual multipliers; returns true on success
template <class R>
bool SoPlexBase<R>::getDualViolation(R& maxviol, R& sumviol)
{
  if( !isDualFeasible() || !hasBasis() )
    return false;

  _syncRealSolution();
  VectorBase<R>& dual = _solReal._dual;
  assert(dual.dim() == numRows());

  maxviol = 0.0;
  sumviol = 0.0;

  for( int r = numRows() - 1; r >= 0; r-- )
    {
      typename SPxSolverBase<R>::VarStatus rowStatus = basisRowStatus(r);

      if( intParam(SoPlexBase<R>::OBJSENSE) == OBJSENSE_MINIMIZE )
        {
          if( rowStatus != SPxSolverBase<R>::ON_UPPER && rowStatus != SPxSolverBase<R>::FIXED && dual[r] < 0.0 )
            {
              sumviol += -dual[r];
              if( dual[r] < -maxviol )
                maxviol = -dual[r];
            }
          if( rowStatus != SPxSolverBase<R>::ON_LOWER && rowStatus != SPxSolverBase<R>::FIXED && dual[r] > 0.0 )
            {
              sumviol += dual[r];
              if( dual[r] > maxviol )
                maxviol = dual[r];
            }
        }
      else
        {
          if( rowStatus != SPxSolverBase<R>::ON_UPPER && rowStatus != SPxSolverBase<R>::FIXED && dual[r] > 0.0 )
            {
              sumviol += dual[r];
              if( dual[r] > maxviol )
                maxviol = dual[r];
            }
          if( rowStatus != SPxSolverBase<R>::ON_LOWER && rowStatus != SPxSolverBase<R>::FIXED && dual[r] < 0.0 )
            {
              sumviol += -dual[r];
              if( dual[r] < -maxviol )
                maxviol = -dual[r];
            }
        }
    }

  return true;
}

/// gets violation of reduced costs; returns true on success
template <class R>
bool SoPlexBase<R>::getRedCostViolation(R& maxviol, R& sumviol)
{
  if( !isDualFeasible() || !hasBasis() )
    return false;

  _syncRealSolution();
  VectorBase<R>& redcost = _solReal._redCost;
  assert(redcost.dim() == numCols());

  maxviol = 0.0;
  sumviol = 0.0;

  for( int c = numCols() - 1; c >= 0; c-- )
    {
      typename SPxSolverBase<R>::VarStatus colStatus = basisColStatus(c);

      if( intParam(SoPlexBase<R>::OBJSENSE) == OBJSENSE_MINIMIZE )
        {
          if( colStatus != SPxSolverBase<R>::ON_UPPER && colStatus != SPxSolverBase<R>::FIXED && redcost[c] < 0.0 )
            {
              sumviol += -redcost[c];
              if( redcost[c] < -maxviol )
                maxviol = -redcost[c];
            }
          if( colStatus != SPxSolverBase<R>::ON_LOWER && colStatus != SPxSolverBase<R>::FIXED && redcost[c] > 0.0 )
            {
              sumviol += redcost[c];
              if( redcost[c] > maxviol )
                maxviol = redcost[c];
            }
        }
      else
        {
          if( colStatus != SPxSolverBase<R>::ON_UPPER && colStatus != SPxSolverBase<R>::FIXED && redcost[c] > 0.0 )
            {
              sumviol += redcost[c];
              if( redcost[c] > maxviol )
                maxviol = redcost[c];
            }
          if( colStatus != SPxSolverBase<R>::ON_LOWER && colStatus != SPxSolverBase<R>::FIXED && redcost[c] < 0.0 )
            {
              sumviol += -redcost[c];
              if( redcost[c] < -maxviol )
                maxviol = -redcost[c];
            }
        }
    }

  return true;
}


/// assignment operator
template <class R>
SoPlexBase<R>& SoPlexBase<R>::operator=(const SoPlexBase<R>& rhs)
{
  if( this != &rhs )
    {
      // copy message handler
      spxout = rhs.spxout;

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
      _scalerGeoequi = rhs._scalerGeoequi;
      _scalerLeastsq = rhs._scalerLeastsq;
      _starterWeight = rhs._starterWeight;
      _starterSum = rhs._starterSum;
      _starterVector = rhs._starterVector;
      _pricerAuto = rhs._pricerAuto;
      _pricerDantzig = rhs._pricerDantzig;
      _pricerParMult = rhs._pricerParMult;
      _pricerDevex = rhs._pricerDevex;
      _pricerQuickSteep = rhs._pricerQuickSteep;
      _pricerSteep = rhs._pricerSteep;
      _ratiotesterTextbook = rhs._ratiotesterTextbook;
      _ratiotesterHarris = rhs._ratiotesterHarris;
      _ratiotesterFast = rhs._ratiotesterFast;
      _ratiotesterBoundFlipping = rhs._ratiotesterBoundFlipping;

      // copy solution data
      _status = rhs._status;
      _lastSolveMode = rhs._lastSolveMode;
      _basisStatusRows = rhs._basisStatusRows;
      _basisStatusCols = rhs._basisStatusCols;

      if( rhs._hasSolReal )
        _solReal = rhs._solReal;

      if( rhs._hasSolRational )
        _solRational = rhs._solRational;

      // set message handlers in members
      _solver.setOutstream(spxout);
      _scalerUniequi.setOutstream(spxout);
      _scalerBiequi.setOutstream(spxout);
      _scalerGeo1.setOutstream(spxout);
      _scalerGeo8.setOutstream(spxout);
      _scalerGeoequi.setOutstream(spxout);
      _scalerLeastsq.setOutstream(spxout);

      // transfer the lu solver
      _solver.setBasisSolver(&_slufactor);

      // initialize pointers for simplifier, scaler, and starter
      setIntParam(SoPlexBase<R>::SIMPLIFIER, intParam(SoPlexBase<R>::SIMPLIFIER), true);
      setIntParam(SoPlexBase<R>::SCALER, intParam(SoPlexBase<R>::SCALER), true);
      setIntParam(SoPlexBase<R>::STARTER, intParam(SoPlexBase<R>::STARTER), true);

      // copy real LP if different from the LP in the solver
      if( rhs._realLP != &(rhs._solver) )
        {
          _realLP = 0;
          spx_alloc(_realLP);
          _realLP = new (_realLP) SPxLPBase<R>(*(rhs._realLP));
        }
      else
        _realLP = &_solver;

      // copy rational LP
      if( rhs._rationalLP == 0 )
        {
          assert(intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL);
          _rationalLP = 0;
        }
      else
        {
          assert(intParam(SoPlexBase<R>::SYNCMODE) != SYNCMODE_ONLYREAL);
          _rationalLP = 0;
          spx_alloc(_rationalLP);
          _rationalLP = new (_rationalLP) SPxLPRational(*rhs._rationalLP);
        }

      // copy rational factorization
      if( rhs._rationalLUSolver.status() == SLinSolverRational::OK )
        _rationalLUSolver = rhs._rationalLUSolver;

      // copy boolean flags
      _isRealLPLoaded = rhs._isRealLPLoaded;
      _isRealLPScaled = rhs._isRealLPScaled;
      _hasSolReal = rhs._hasSolReal;
      _hasSolRational = rhs._hasSolRational;
      _hasBasis = rhs._hasBasis;
      _applyPolishing = rhs._applyPolishing;

      // rational constants do not need to be assigned
      _rationalPosone = 1;
      _rationalNegone = -1;
      _rationalZero = 0;
    }

  assert(_isConsistent());

  return *this;
}

/// returns smallest non-zero element in absolute value
template <class R>
R SoPlexBase<R>::minAbsNonzeroReal() const
{
  assert(_realLP != 0);
  return _realLP->minAbsNzo();
}


/// returns biggest non-zero element in absolute value
template <class R>
R SoPlexBase<R>::maxAbsNonzeroReal() const
{
  assert(_realLP != 0);
  return _realLP->maxAbsNzo();
}


/// returns (unscaled) coefficient
template <class R>
R SoPlexBase<R>::coefReal(int row, int col) const
{
  if( _realLP->isScaled() )
    {
      assert(_scaler);
      return _scaler->getCoefUnscaled(*_realLP, row, col);
    }
  else
    return colVectorRealInternal(col)[row];
}

/// returns vector of row \p i, ignoring scaling
template <class R>
const SVectorBase<R>& SoPlexBase<R>::rowVectorRealInternal(int i) const
{
  assert(_realLP != 0);
  return _realLP->rowVector(i);
}

/// gets vector of row \p i
template <class R>
void SoPlexBase<R>::getRowVectorReal(int i, DSVectorBase<R>& row) const
{
  assert(_realLP);

  if( _realLP->isScaled() )
    {
      assert(_scaler);
      row.setMax(_realLP->rowVector(i).size());
      _scaler->getRowUnscaled(*_realLP, i, row);
    }
  else
    row = _realLP->rowVector(i);
}


/// returns right-hand side vector, ignoring scaling
template <class R>
const VectorBase<R>& SoPlexBase<R>::rhsRealInternal() const
{
  assert(_realLP != 0);
  return _realLP->rhs();
}



/// gets right-hand side vector
template <class R>
void SoPlexBase<R>::getRhsReal(DVectorBase<R>& rhs) const
{
  assert(_realLP);
  _realLP->getRhsUnscaled(rhs);
}



/// returns right-hand side of row \p i
template <class R>
R SoPlexBase<R>::rhsReal(int i) const
{
  assert(_realLP != 0);
  return _realLP->rhsUnscaled(i);
}

/// returns left-hand side vector, ignoring scaling
template <class R>
const VectorBase<R>& SoPlexBase<R>::lhsRealInternal() const
{
  assert(_realLP != 0);
  return _realLP->lhs();
}

/// gets left-hand side vector
template <class R>
void SoPlexBase<R>::getLhsReal(DVectorBase<R>& lhs) const
{
  assert(_realLP);
  _realLP->getLhsUnscaled(lhs);
}

/// returns left-hand side of row \p i
template <class R>
R SoPlexBase<R>::lhsReal(int i) const
{
  assert(_realLP != 0);
  return _realLP->lhsUnscaled(i);
}


/// returns inequality type of row \p i
template <class R>
LPRowReal::Type SoPlexBase<R>::rowTypeReal(int i) const
{
  assert(_realLP != 0);
  return _realLP->rowType(i);
}

/// returns vector of col \p i, ignoring scaling
template <class R>
const SVectorBase<R>& SoPlexBase<R>::colVectorRealInternal(int i) const
{
  assert(_realLP != 0);
  return _realLP->colVector(i);
}

/// gets vector of col \p i
template <class R>
void SoPlexBase<R>::getColVectorReal(int i, DSVectorBase<R>& col) const
{
  assert(_realLP);
  _realLP->getColVectorUnscaled(i, col);
}


/// returns upper bound vector
template <class R>
const VectorBase<R>& SoPlexBase<R>::upperRealInternal() const
{
  assert(_realLP != 0);
  return _realLP->upper();
}


/// returns upper bound of column \p i
template <class R>
R SoPlexBase<R>::upperReal(int i) const
{
  assert(_realLP != 0);
  return _realLP->upperUnscaled(i);
}


/// gets upper bound vector
template <class R>
void SoPlexBase<R>::getUpperReal(DVectorBase<R>& upper) const
{
  assert(_realLP != 0);
  return _realLP->getUpperUnscaled(upper);
}


/// returns lower bound vector
template <class R>
const VectorBase<R>& SoPlexBase<R>::lowerRealInternal() const
{
  assert(_realLP != 0);
  return _realLP->lower();
}



/// returns lower bound of column \p i
template <class R>
R SoPlexBase<R>::lowerReal(int i) const
{
  assert(_realLP != 0);
  return _realLP->lowerUnscaled(i);
}


/// gets lower bound vector
template <class R>
void SoPlexBase<R>::getLowerReal(DVectorBase<R>& lower) const
{
  assert(_realLP != 0);
  return _realLP->getLowerUnscaled(lower);
}


/// gets objective function vector
template <class R>
void SoPlexBase<R>::getObjReal(VectorBase<R>& obj) const
{
  assert(_realLP != 0);
  _realLP->getObjUnscaled(obj);
}


/// returns objective value of column \p i
template <class R>
R SoPlexBase<R>::objReal(int i) const
{
  assert(_realLP != 0);
  return _realLP->objUnscaled(i);
}


/// returns objective function vector after transformation to a maximization problem; since this is how it is stored
/// internally, this is generally faster
template <class R>
const VectorBase<R>& SoPlexBase<R>::maxObjRealInternal() const
{
  assert(_realLP != 0);
  return _realLP->maxObj();
}


/// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
/// stored internally, this is generally faster
template <class R>
R SoPlexBase<R>::maxObjReal(int i) const
{
  assert(_realLP != 0);
  return _realLP->maxObjUnscaled(i);
}


/// gets number of available dual norms
template <class R>
void SoPlexBase<R>::getNdualNorms(int& nnormsRow, int& nnormsCol) const
{
  _solver.getNdualNorms(nnormsRow, nnormsCol);
}


/// gets steepest edge norms and returns false if they are not available
template <class R>
bool SoPlexBase<R>::getDualNorms(int& nnormsRow, int& nnormsCol, R* norms) const
{
  return _solver.getDualNorms(nnormsRow, nnormsCol, norms);
}


/// sets steepest edge norms and returns false if that's not possible
template <class R>
bool SoPlexBase<R>::setDualNorms(int nnormsRow, int nnormsCol, R* norms)
{
  return _solver.setDualNorms(nnormsRow, nnormsCol, norms);
}


/// pass integrality information about the variables to the solver
template <class R>
void SoPlexBase<R>::setIntegralityInformation( int ncols, int* intInfo)
{
  assert(ncols == _solver.nCols() || (ncols == 0 && intInfo == NULL));
  _solver.setIntegralityInformation(ncols, intInfo);
}


template <class R>
int SoPlexBase<R>::numRowsRational() const
{
  assert(_rationalLP != 0);
  return _rationalLP->nRows();
}

/// returns number of columns
template <class R>
int SoPlexBase<R>::numColsRational() const
{
  assert(_rationalLP != 0);
  return _rationalLP->nCols();
}

/// returns number of nonzeros
template <class R>
int SoPlexBase<R>::numNonzerosRational() const
{
  assert(_rationalLP != 0);
  return _rationalLP->nNzos();
}

/// returns smallest non-zero element in absolute value
template <class R>
Rational SoPlexBase<R>::minAbsNonzeroRational() const
{
  assert(_rationalLP != 0);
  return _rationalLP->minAbsNzo();
}

/// returns biggest non-zero element in absolute value
template <class R>
Rational SoPlexBase<R>::maxAbsNonzeroRational() const
{
  assert(_rationalLP != 0);
  return _rationalLP->maxAbsNzo();
}

/// gets row \p i
template <class R>
void SoPlexBase<R>::getRowRational(int i, LPRowRational& lprow) const
{
  assert(_rationalLP != 0);
  _rationalLP->getRow(i, lprow);
}


/// gets rows \p start, ..., \p end.
template <class R>
void SoPlexBase<R>::getRowsRational(int start, int end, LPRowSetRational& lprowset) const
{
  assert(_rationalLP != 0);
  _rationalLP->getRows(start, end, lprowset);
}



/// returns vector of row \p i
template <class R>
const SVectorRational& SoPlexBase<R>::rowVectorRational(int i) const
{
  assert(_rationalLP != 0);
  return _rationalLP->rowVector(i);
}

/// returns right-hand side vector
template <class R>
const VectorRational& SoPlexBase<R>::rhsRational() const
{
  assert(_rationalLP != 0);
  return _rationalLP->rhs();
}



/// returns right-hand side of row \p i
template <class R>
const Rational& SoPlexBase<R>::rhsRational(int i) const
{
  assert(_rationalLP != 0);
  return _rationalLP->rhs(i);
}


/// returns left-hand side vector
template <class R>
const VectorRational& SoPlexBase<R>::lhsRational() const
{
  assert(_rationalLP != 0);
  return _rationalLP->lhs();
}



/// returns left-hand side of row \p i
template <class R>
const Rational& SoPlexBase<R>::lhsRational(int i) const
{
  assert(_rationalLP != 0);
  return _rationalLP->lhs(i);
}



/// returns inequality type of row \p i
template <class R>
LPRowRational::Type SoPlexBase<R>::rowTypeRational(int i) const
{
  assert(_rationalLP != 0);
  return _rationalLP->rowType(i);
}



/// gets column \p i
template <class R>
void SoPlexBase<R>::getColRational(int i, LPColRational& lpcol) const
{
  assert(_rationalLP != 0);
  return _rationalLP->getCol(i, lpcol);
}



/// gets columns \p start, ..., \p end
template <class R>
void SoPlexBase<R>::getColsRational(int start, int end, LPColSetRational& lpcolset) const
{
  assert(_rationalLP != 0);
  return _rationalLP->getCols(start, end, lpcolset);
}


/// returns vector of column \p i
template <class R>
const SVectorRational& SoPlexBase<R>::colVectorRational(int i) const
{
  assert(_rationalLP != 0);
  return _rationalLP->colVector(i);
}



/// returns upper bound vector
template <class R>
const VectorRational& SoPlexBase<R>::upperRational() const
{
  assert(_rationalLP != 0);
  return _rationalLP->upper();
}



/// returns upper bound of column \p i
template <class R>
const Rational& SoPlexBase<R>::upperRational(int i) const
{
  assert(_rationalLP != 0);
  return _rationalLP->upper(i);
}



/// returns lower bound vector
template <class R>
const VectorRational& SoPlexBase<R>::lowerRational() const
{
  assert(_rationalLP != 0);
  return _rationalLP->lower();
}

/// returns lower bound of column \p i
template <class R>
const Rational& SoPlexBase<R>::lowerRational(int i) const
{
  assert(_rationalLP != 0);
  return _rationalLP->lower(i);
}



/// gets objective function vector
template <class R>
void SoPlexBase<R>::getObjRational(VectorRational& obj) const
{
  assert(_rationalLP != 0);
  _rationalLP->getObj(obj);
}



/// gets objective value of column \p i
template <class R>
void SoPlexBase<R>::getObjRational(int i, Rational& obj) const
{
  obj = maxObjRational(i);
  if( intParam(SoPlexBase<R>::OBJSENSE) == SoPlexBase<R>::OBJSENSE_MINIMIZE )
    obj *= -1;
}



/// returns objective value of column \p i
template <class R>
Rational SoPlexBase<R>::objRational(int i) const
{
  assert(_rationalLP != 0);
  return _rationalLP->obj(i);
}



/// returns objective function vector after transformation to a maximization problem; since this is how it is stored
/// internally, this is generally faster
template <class R>
const VectorRational& SoPlexBase<R>::maxObjRational() const
{
  assert(_rationalLP != 0);
  return _rationalLP->maxObj();
}



/// returns objective value of column \p i after transformation to a maximization problem; since this is how it is
/// stored internally, this is generally faster
template <class R>
const Rational& SoPlexBase<R>::maxObjRational(int i) const
{
  assert(_rationalLP != 0);
  return _rationalLP->maxObj(i);
}



/// adds a single row
template <class R>
void SoPlexBase<R>::addRowReal(const LPRowReal& lprow)
{
  assert(_realLP != 0);

  _addRowReal(lprow);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->addRow(lprow);
      _completeRangeTypesRational();
    }

  _invalidateSolution();
}



/// adds multiple rows
template <class R>
void SoPlexBase<R>::addRowsReal(const LPRowSetReal& lprowset)
{
  assert(_realLP != 0);

  _addRowsReal(lprowset);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->addRows(lprowset);
      _completeRangeTypesRational();
    }

  _invalidateSolution();
}



/// adds a single column
template <class R>
void SoPlexBase<R>::addColReal(const LPColReal& lpcol)
{
  assert(_realLP != 0);

  _addColReal(lpcol);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->addCol(lpcol);
      _completeRangeTypesRational();
    }

  _invalidateSolution();
}



/// adds multiple columns
template <class R>
void SoPlexBase<R>::addColsReal(const LPColSetReal& lpcolset)
{
  assert(_realLP != 0);

  _addColsReal(lpcolset);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->addCols(lpcolset);
      _completeRangeTypesRational();
    }

  _invalidateSolution();
}



/// replaces row \p i with \p lprow
template <class R>
void SoPlexBase<R>::changeRowReal(int i, const LPRowReal& lprow)
{
  assert(_realLP != 0);

  _changeRowReal(i, lprow);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeRow(i, lprow);
      _rowTypes[i] = _rangeTypeReal(lprow.lhs(), lprow.rhs());
      _completeRangeTypesRational();
    }

  _invalidateSolution();
}
#ifdef SOPLEX_WITH_GMP
/// adds a single column
template <class R>
void SoPlexBase<R>::addColRational(const mpq_t* obj, const mpq_t* lower, const mpq_t* colValues, const int* colIndices, const int colSize, const mpq_t* upper)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->addCol(obj, lower, colValues, colIndices, colSize, upper);
  int i = numColsRational() - 1;
  _completeRangeTypesRational();

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _addColReal(Real(maxObjRational(i)) * (intParam(SoPlexBase<R>::OBJSENSE) == SoPlexBase<R>::OBJSENSE_MAXIMIZE ? 1.0 : -1.0),
                R(lowerRational(i)), DSVectorBase<R>(_rationalLP->colVector(i)), R(upperRational(i)));

  _invalidateSolution();
}



/// adds a set of columns
template <class R>
void SoPlexBase<R>::addColsRational(const mpq_t* obj, const mpq_t* lower, const mpq_t* colValues, const int* colIndices, const int* colStarts, const int* colLengths, const int numCols, const int numValues, const mpq_t* upper)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->addCols(obj, lower, colValues, colIndices, colStarts, colLengths, numCols, numValues, upper);
  _completeRangeTypesRational();

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      LPColSetReal lpcolset;
      for( int i = numColsRational() - numCols; i < numColsRational(); i++ )
        lpcolset.add(Real(maxObjRational(i)) * (intParam(SoPlexBase<R>::OBJSENSE) == SoPlexBase<R>::OBJSENSE_MAXIMIZE ? 1.0 : -1.0),
                     R(lowerRational(i)), DSVectorBase<R>(_rationalLP->colVector(i)), R(upperRational(i)));
      _addColsReal(lpcolset);
    }

  _invalidateSolution();
}
#endif


/// changes left-hand side vector for constraints to \p lhs
template <class R>
void SoPlexBase<R>::changeLhsReal(const VectorBase<R>& lhs)
{
  assert(_realLP != 0);

  _changeLhsReal(lhs);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeLhs(DVectorRational(lhs));
      for( int i = 0; i < numRowsRational(); i++ )
        _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
    }

  _invalidateSolution();
}



/// changes left-hand side of row \p i to \p lhs
template <class R>
void SoPlexBase<R>::changeLhsReal(int i, const R& lhs)
{
  assert(_realLP != 0);

  _changeLhsReal(i, lhs);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeLhs(i, lhs);
      _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
    }

  _invalidateSolution();
}



/// changes right-hand side vector to \p rhs
template <class R>
void SoPlexBase<R>::changeRhsReal(const VectorBase<R>& rhs)
{
  assert(_realLP != 0);

  _changeRhsReal(rhs);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeRhs(DVectorRational(rhs));
      for( int i = 0; i < numRowsRational(); i++ )
        _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
    }

  _invalidateSolution();
}



/// changes right-hand side of row \p i to \p rhs
template <class R>
void SoPlexBase<R>::changeRhsReal(int i, const R& rhs)
{
  assert(_realLP != 0);

  _changeRhsReal(i, rhs);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeRhs(i, rhs);
      _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
    }

  _invalidateSolution();
}


/// changes left- and right-hand side vectors
template <class R>
void SoPlexBase<R>::changeRangeReal(const VectorBase<R>& lhs, const VectorBase<R>& rhs)
{
  assert(_realLP != 0);

  _changeRangeReal(lhs, rhs);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeRange(DVectorRational(lhs), DVectorRational(rhs));
      for( int i = 0; i < numRowsRational(); i++ )
        _rowTypes[i] = _rangeTypeReal(lhs[i], rhs[i]);
    }

  _invalidateSolution();
}



/// changes left- and right-hand side of row \p i
template <class R>
void SoPlexBase<R>::changeRangeReal(int i, const R& lhs, const R& rhs)
{
  assert(_realLP != 0);

  _changeRangeReal(i,lhs, rhs);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeRange(i, lhs, rhs);
      _rowTypes[i] = _rangeTypeReal(lhs, rhs);
    }

  _invalidateSolution();
}



/// replaces column \p i with \p lpcol
template <class R>
void SoPlexBase<R>::changeColReal(int i, const LPColReal& lpcol)
{
  assert(_realLP != 0);

  _changeColReal(i, lpcol);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeCol(i, lpcol);
      _colTypes[i] = _rangeTypeReal(lpcol.lower(), lpcol.upper());
      _completeRangeTypesRational();
    }

  _invalidateSolution();
}



/// changes vector of lower bounds to \p lower
template <class R>
void SoPlexBase<R>::changeLowerReal(const VectorBase<R>& lower)
{
  assert(_realLP != 0);

  _changeLowerReal(lower);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeLower(DVectorRational(lower));
      for( int i = 0; i < numColsRational(); i++ )
        _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
    }


  _invalidateSolution();
}



/// changes lower bound of column i to \p lower
template <class R>
void SoPlexBase<R>::changeLowerReal(int i, const R& lower)
{
  assert(_realLP != 0);

  _changeLowerReal(i, lower);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeLower(i, lower);
      _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
    }

  _invalidateSolution();
}



/// changes vector of upper bounds to \p upper
template <class R>
void SoPlexBase<R>::changeUpperReal(const VectorBase<R>& upper)
{
  assert(_realLP != 0);

  _changeUpperReal(upper);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeUpper(DVectorRational(upper));
      for( int i = 0; i < numColsRational(); i++ )
        _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
    }

  _invalidateSolution();
}



/// changes \p i 'th upper bound to \p upper
template <class R>
void SoPlexBase<R>::changeUpperReal(int i, const R& upper)
{
  assert(_realLP != 0);

  _changeUpperReal(i, upper);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeUpper(i, upper);
      _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));
    }

  _invalidateSolution();
}



/// changes vectors of column bounds to \p lower and \p upper
template <class R>
void SoPlexBase<R>::changeBoundsReal(const VectorBase<R>& lower, const VectorBase<R>& upper)
{
  assert(_realLP != 0);

  _changeBoundsReal(lower, upper);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeBounds(DVectorRational(lower), DVectorRational(upper));
      for( int i = 0; i < numColsRational(); i++ )
        _colTypes[i] = _rangeTypeReal(lower[i], upper[i]);
    }

  _invalidateSolution();
}



/// changes bounds of column \p i to \p lower and \p upper
template <class R>
void SoPlexBase<R>::changeBoundsReal(int i, const R& lower, const R& upper)
{
  assert(_realLP != 0);

  _changeBoundsReal(i, lower, upper);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->changeBounds(i, lower, upper);
      _colTypes[i] = _rangeTypeReal(lower, upper);
    }
  _invalidateSolution();
}


/// changes objective function vector to \p obj
template <class R>
void SoPlexBase<R>::changeObjReal(const VectorBase<R>& obj)
{
  assert(_realLP != 0);

  bool scale = _realLP->isScaled();
  _realLP->changeObj(obj, scale);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _rationalLP->changeObj(DVectorRational(obj));

  _invalidateSolution();
}



/// changes objective coefficient of column i to \p obj
template <class R>
void SoPlexBase<R>::changeObjReal(int i, const R& obj)
{
  assert(_realLP != 0);

  bool scale = _realLP->isScaled();
  _realLP->changeObj(i, obj, scale);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _rationalLP->changeObj(i, obj);

  _invalidateSolution();
}



/// changes matrix entry in row \p i and column \p j to \p val
template <class R>
void SoPlexBase<R>::changeElementReal(int i, int j, const R& val)
{
  assert(_realLP != 0);

  _changeElementReal(i, j, val);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _rationalLP->changeElement(i, j, val);

  _invalidateSolution();
}



/// removes row \p i
template <class R>
void SoPlexBase<R>::removeRowReal(int i)
{
  assert(_realLP != 0);

  _removeRowReal(i);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->removeRow(i);
      // only swap elements if not the last one was removed
      if( i < _rationalLP->nRows() )
        {
          _rowTypes[i] = _rowTypes[_rationalLP->nRows()];
          assert(_rowTypes[i] == _rangeTypeRational(lhsRational(i), rhsRational(i)));
        }
      _rowTypes.reSize(_rationalLP->nRows());
    }

  _invalidateSolution();
}



/// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
/// new index where row \p i has been moved to; note that \p perm must point to an array of size at least
/// #numRows()
template <class R>
void SoPlexBase<R>::removeRowsReal(int perm[])
{
  assert(_realLP != 0);

  const int oldsize = numRows();
  _removeRowsReal(perm);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->removeRows(perm);
      for( int i = 0; i < oldsize; i++ )
        {
          if( perm[i] >= 0 )
            _rowTypes[perm[i]] = _rowTypes[i];
        }
      _rowTypes.reSize(_rationalLP->nRows());
      for( int i = 0; i < numRowsRational(); i++ )
        {
          assert(_rowTypes[i] == _rangeTypeRational(lhsRational(i), rhsRational(i)));
        }
    }

  _invalidateSolution();
}



/// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRows() may be passed
/// as buffer memory
template <class R>
void SoPlexBase<R>::removeRowsReal(int idx[], int n, int perm[])
{
  if( perm == 0 )
    {
      DataArray< int > p(numRows());
      _idxToPerm(idx, n, p.get_ptr(), numRows());
      SoPlexBase<R>::removeRowsReal(p.get_ptr());
    }
  else
    {
      _idxToPerm(idx, n, perm, numRows());
      SoPlexBase<R>::removeRowsReal(perm);
    }
}



/// removes rows \p start to \p end including both; an array \p perm of size #numRows() may be passed as buffer
/// memory
template <class R>
void SoPlexBase<R>::removeRowRangeReal(int start, int end, int perm[])
{
  if( perm == 0 )
    {
      DataArray< int > p(numRows());
      _rangeToPerm(start, end, p.get_ptr(), numRows());
      SoPlexBase<R>::removeRowsReal(p.get_ptr());
    }
  else
    {
      _rangeToPerm(start, end, perm, numRows());
      SoPlexBase<R>::removeRowsReal(perm);
    }
}



/// removes column i
template <class R>
void SoPlexBase<R>::removeColReal(int i)
{
  assert(_realLP != 0);

  _removeColReal(i);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->removeCol(i);
      // only swap elements if not the last one was removed
      if( i < _rationalLP->nCols() )
        {
          _colTypes[i] = _colTypes[_rationalLP->nCols()];
          assert(_colTypes[i] == _rangeTypeRational(lowerRational(i), upperRational(i)));
        }
      _colTypes.reSize(_rationalLP->nCols());
    }

  _invalidateSolution();
}



/// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
/// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
/// #numColsReal()
template <class R>
void SoPlexBase<R>::removeColsReal(int perm[])
{
  assert(_realLP != 0);

  const int oldsize = numCols();
  _removeColsReal(perm);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->removeCols(perm);
      for( int i = 0; i < oldsize; i++ )
        {
          if( perm[i] >= 0 )
            _colTypes[perm[i]] = _colTypes[i];
        }
      _colTypes.reSize(_rationalLP->nCols());
      for( int i = 0; i < numColsRational(); i++ )
        {
          assert(_colTypes[i] == _rangeTypeRational(lowerRational(i), upperRational(i)));
        }
    }

  _invalidateSolution();
}



/// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsReal() may be
/// passed as buffer memory
template <class R>
void SoPlexBase<R>::removeColsReal(int idx[], int n, int perm[])
{
  if( perm == 0 )
    {
      DataArray< int > p(numCols());
      _idxToPerm(idx, n, p.get_ptr(), numCols());
      SoPlexBase<R>::removeColsReal(p.get_ptr());
    }
  else
    {
      _idxToPerm(idx, n, perm, numCols());
      SoPlexBase<R>::removeColsReal(perm);
    }
}



/// removes columns \p start to \p end including both; an array \p perm of size #numCols() may be passed as
/// buffer memory
template <class R>
void SoPlexBase<R>::removeColRangeReal(int start, int end, int perm[])
{
  if( perm == 0 )
    {
      DataArray< int > p(numCols());
      _rangeToPerm(start, end, p.get_ptr(), numCols());
      SoPlexBase<R>::removeColsReal(p.get_ptr());
    }
  else
    {
      _rangeToPerm(start, end, perm, numCols());
      SoPlexBase<R>::removeColsReal(perm);
    }
}



/// clears the LP
template <class R>
void SoPlexBase<R>::clearLPReal()
{
  assert(_realLP != 0);

  _realLP->clear();
  _hasBasis = false;
  _rationalLUSolver.clear();

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      _rationalLP->clear();
      _rowTypes.clear();
      _colTypes.clear();
    }

  _invalidateSolution();
}



/// synchronizes R LP with rational LP, i.e., copies (rounded) rational LP into R LP, if sync mode is manual
template <class R>
void SoPlexBase<R>::syncLPReal()
{
  assert(_isConsistent());

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_MANUAL )
    _syncLPReal();
}



/// adds a single row
template <class R>
void SoPlexBase<R>::addRowRational(const LPRowRational& lprow)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->addRow(lprow);
  _completeRangeTypesRational();

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _addRowReal(lprow);

  _invalidateSolution();
}



#ifdef SOPLEX_WITH_GMP
/// adds a single row
template <class R>
void SoPlexBase<R>::addRowRational(const mpq_t* lhs, const mpq_t* rowValues, const int* rowIndices, const int rowSize, const mpq_t* rhs)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->addRow(lhs, rowValues, rowIndices, rowSize, rhs);
  _completeRangeTypesRational();

  int i = numRowsRational() - 1;
  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _addRowReal(Real(lhsRational(i)), DSVectorBase<R>(_rationalLP->rowVector(i)), R(rhsRational(i)));

  _invalidateSolution();
}



/// adds a set of rows
template <class R>
void SoPlexBase<R>::addRowsRational(const mpq_t* lhs, const mpq_t* rowValues, const int* rowIndices, const int* rowStarts, const int* rowLengths, const int numRows, const int numValues, const mpq_t* rhs)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->addRows(lhs, rowValues, rowIndices, rowStarts, rowLengths, numRows, numValues, rhs);
  _completeRangeTypesRational();

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    {
      LPRowSetReal lprowset;
      for( int i = numRowsRational() - numRows; i < numRowsRational(); i++ )
        lprowset.add(Real(lhsRational(i)), DSVectorBase<R>(_rationalLP->rowVector(i)), R(rhsRational(i)));
      _addRowsReal(lprowset);
    }

  _invalidateSolution();
}
#endif



/// adds multiple rows
template <class R>
void SoPlexBase<R>::addRowsRational(const LPRowSetRational& lprowset)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->addRows(lprowset);
  _completeRangeTypesRational();

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _addRowsReal(lprowset);

  _invalidateSolution();
}


/// adds a single column
template <class R>
void SoPlexBase<R>::addColRational(const LPColRational& lpcol)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->addCol(lpcol);
  _completeRangeTypesRational();

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _addColReal(lpcol);

  _invalidateSolution();
}





/// adds multiple columns
template <class R>
void SoPlexBase<R>::addColsRational(const LPColSetRational& lpcolset)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->addCols(lpcolset);
  _completeRangeTypesRational();

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _addColsReal(lpcolset);

  _invalidateSolution();
}



/// replaces row \p i with \p lprow
template <class R>
void SoPlexBase<R>::changeRowRational(int i, const LPRowRational& lprow)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeRow(i, lprow);
  _rowTypes[i] = _rangeTypeRational(lprow.lhs(), lprow.rhs());
  _completeRangeTypesRational();

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeRowReal(i, lprow);

  _invalidateSolution();
}



/// changes left-hand side vector for constraints to \p lhs
template <class R>
void SoPlexBase<R>::changeLhsRational(const VectorRational& lhs)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeLhs(lhs);
  for( int i = 0; i < numRowsRational(); i++ )
    _rowTypes[i] = _rangeTypeRational(lhs[i], _rationalLP->rhs(i));

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeLhsReal(DVectorBase<R>(lhs));

  _invalidateSolution();
}



/// changes left-hand side of row \p i to \p lhs
template <class R>
void SoPlexBase<R>::changeLhsRational(int i, const Rational& lhs)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeLhs(i, lhs);
  _rowTypes[i] = _rangeTypeRational(lhs, _rationalLP->rhs(i));

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeLhsReal(i, R(lhs));

  _invalidateSolution();
}



#ifdef SOPLEX_WITH_GMP
/// changes left-hand side of row \p i to \p lhs
template <class R>
void SoPlexBase<R>::changeLhsRational(int i, const mpq_t* lhs)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeLhs(i, lhs);
  _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeLhsReal(i, R(lhsRational(i)));

  _invalidateSolution();
}
#endif



/// changes right-hand side vector to \p rhs
template <class R>
void SoPlexBase<R>::changeRhsRational(const VectorRational& rhs)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeRhs(rhs);
  for( int i = 0; i < numRowsRational(); i++ )
    _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), rhs[i]);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeRhsReal(DVectorBase<R>(rhs));

  _invalidateSolution();
}



#ifdef SOPLEX_WITH_GMP
/// changes right-hand side vector to \p rhs
template <class R>
void SoPlexBase<R>::changeRhsRational(const mpq_t* rhs, int rhsSize)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  for( int i = 0; i < rhsSize; i++ )
    {
      _rationalLP->changeRhs(i, rhs[i]);
      _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));
    }

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeRhsReal(DVectorBase<R>(rhsRational()));

  _invalidateSolution();
}
#endif



/// changes right-hand side of row \p i to \p rhs
template <class R>
void SoPlexBase<R>::changeRhsRational(int i, const Rational& rhs)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeRhs(i, rhs);
  _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), rhs);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeRhsReal(i, R(rhs));

  _invalidateSolution();
}



/// changes left- and right-hand side vectors
template <class R>
void SoPlexBase<R>::changeRangeRational(const VectorRational& lhs, const VectorRational& rhs)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeRange(lhs, rhs);
  for( int i = 0; i < numRowsRational(); i++ )
    _rowTypes[i] = _rangeTypeRational(lhs[i], rhs[i]);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeRangeReal(DVectorBase<R>(lhs), DVectorBase<R>(rhs));

  _invalidateSolution();
}



/// changes left- and right-hand side of row \p i
template <class R>
void SoPlexBase<R>::changeRangeRational(int i, const Rational& lhs, const Rational& rhs)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeRange(i, lhs, rhs);
  _rowTypes[i] = _rangeTypeRational(lhs, rhs);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeRangeReal(i, R(lhs), R(rhs));

  _invalidateSolution();
}



#ifdef SOPLEX_WITH_GMP
/// changes left-hand side of row \p i to \p lhs
template <class R>
void SoPlexBase<R>::changeRangeRational(int i, const mpq_t* lhs, const mpq_t* rhs)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeRange(i, lhs, rhs);
  _rowTypes[i] = _rangeTypeRational(_rationalLP->lhs(i), _rationalLP->rhs(i));

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeRangeReal(i, R(lhsRational(i)), R(rhsRational(i)));

  _invalidateSolution();
}
#endif



/// replaces column \p i with \p lpcol
template <class R>
void SoPlexBase<R>::changeColRational(int i, const LPColRational& lpcol)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeCol(i, lpcol);
  _colTypes[i] = _rangeTypeRational(lpcol.lower(), lpcol.upper());
  _completeRangeTypesRational();

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeColReal(i, lpcol);

  _invalidateSolution();
}



/// changes vector of lower bounds to \p lower
template <class R>
void SoPlexBase<R>::changeLowerRational(const VectorRational& lower)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeLower(lower);
  for( int i = 0; i < numColsRational(); i++ )
    _colTypes[i] = _rangeTypeRational(lower[i], _rationalLP->upper(i));

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeLowerReal(DVectorBase<R>(lower));

  _invalidateSolution();
}



/// changes lower bound of column i to \p lower
template <class R>
void SoPlexBase<R>::changeLowerRational(int i, const Rational& lower)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeLower(i, lower);
  _colTypes[i] = _rangeTypeRational(lower, _rationalLP->upper(i));

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeLowerReal(i, R(lower));

  _invalidateSolution();
}



#ifdef SOPLEX_WITH_GMP
/// changes lower bound of column i to \p lower
template <class R>
void SoPlexBase<R>::changeLowerRational(int i, const mpq_t* lower)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeLower(i, lower);
  _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeLowerReal(i, R(lowerRational(i)));

  _invalidateSolution();
}
#endif



/// changes vector of upper bounds to \p upper
template <class R>
void SoPlexBase<R>::changeUpperRational(const VectorRational& upper)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeUpper(upper);
  for( int i = 0; i < numColsRational(); i++ )
    _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), upper[i]);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeUpperReal(DVectorBase<R>(upper));

  _invalidateSolution();
}




/// changes \p i 'th upper bound to \p upper
template <class R>
void SoPlexBase<R>::changeUpperRational(int i, const Rational& upper)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeUpper(i, upper);
  _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), upper);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeUpperReal(i, R(upper));

  _invalidateSolution();
}



#ifdef SOPLEX_WITH_GMP
/// changes upper bound of column i to \p upper
template <class R>
void SoPlexBase<R>::changeUpperRational(int i, const mpq_t* upper)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeUpper(i, upper);
  _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeUpperReal(i, R(upperRational(i)));

  _invalidateSolution();
}
#endif



/// changes vectors of column bounds to \p lower and \p upper
template <class R>
void SoPlexBase<R>::changeBoundsRational(const VectorRational& lower, const VectorRational& upper)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeBounds(lower, upper);
  for( int i = 0; i < numColsRational(); i++ )
    _colTypes[i] = _rangeTypeRational(lower[i], upper[i]);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeBoundsReal(DVectorBase<R>(lower), DVectorBase<R>(upper));

  _invalidateSolution();
}



/// changes bounds of column \p i to \p lower and \p upper
template <class R>
void SoPlexBase<R>::changeBoundsRational(int i, const Rational& lower, const Rational& upper)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeBounds(i, lower, upper);
  _colTypes[i] = _rangeTypeRational(lower, upper);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeBoundsReal(i, R(lower), R(upper));

  _invalidateSolution();
}



#ifdef SOPLEX_WITH_GMP
/// changes bounds of column \p i to \p lower and \p upper
template <class R>
void SoPlexBase<R>::changeBoundsRational(int i, const mpq_t* lower, const mpq_t* upper)
{
  assert(_rationalLP != 0);

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    return;

  _rationalLP->changeBounds(i, lower, upper);
  _colTypes[i] = _rangeTypeRational(_rationalLP->lower(i), _rationalLP->upper(i));

  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
    _changeBoundsReal(i, R(lowerRational(i)), R(upperRational(i)));

  _invalidateSolution();
}
#endif



/// changes objective function vector to \p obj
template <class R>
void SoPlexBase<R>::changeObjRational(const VectorRational& obj)
{
    assert(_rationalLP != 0);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeObj(obj);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
      _realLP->changeObj(DVectorBase<R>(obj));

    _invalidateSolution();
}



/// changes objective coefficient of column i to \p obj
template <class R>
void SoPlexBase<R>::changeObjRational(int i, const Rational& obj)
{
    assert(_rationalLP != 0);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeObj(i, obj);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
      _realLP->changeObj(i, R(obj));

    _invalidateSolution();
}



#ifdef SOPLEX_WITH_GMP
/// changes objective coefficient of column i to \p obj
template <class R>
void SoPlexBase<R>::changeObjRational(int i, const mpq_t* obj)
{
    assert(_rationalLP != 0);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeObj(i, obj);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
      _realLP->changeObj(i, R(objRational(i)));

    _invalidateSolution();
}
#endif



/// changes matrix entry in row \p i and column \p j to \p val
template <class R>
void SoPlexBase<R>::changeElementRational(int i, int j, const Rational& val)
{
    assert(_rationalLP != 0);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeElement(i, j, val);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
      _changeElementReal(i, j, R(val));

    _invalidateSolution();
}


#ifdef SOPLEX_WITH_GMP
/// changes matrix entry in row \p i and column \p j to \p val
template <class R>
void SoPlexBase<R>::changeElementRational(int i, int j, const mpq_t* val)
{
    assert(_rationalLP != 0);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->changeElement(i, j, val);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
      _changeElementReal(i, j, mpq_get_d(*val));

    _invalidateSolution();
}
#endif


/// removes row \p i
template <class R>
void SoPlexBase<R>::removeRowRational(int i)
{
    assert(_rationalLP != 0);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->removeRow(i);
    // only swap elements if not the last one was removed
    if( i < _rationalLP->nRows() )
      {
        _rowTypes[i] = _rowTypes[_rationalLP->nRows()];
        assert(_rowTypes[i] == _rangeTypeRational(lhsRational(i), rhsRational(i)));
      }
    _rowTypes.reSize(_rationalLP->nRows());

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
      _removeRowReal(i);

    _invalidateSolution();
}



/// removes all rows with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the new
/// index where row \p i has been moved to; note that \p perm must point to an array of size at least
/// #numRowsRational()
template <class R>
void SoPlexBase<R>::removeRowsRational(int perm[])
{
    assert(_rationalLP != 0);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    const int oldsize = numRowsRational();
    _rationalLP->removeRows(perm);
    for( int i = 0; i < oldsize; i++ )
      {
        if( perm[i] >= 0 )
          _rowTypes[perm[i]] = _rowTypes[i];
      }
    _rowTypes.reSize(_rationalLP->nRows());
    for( int i = 0; i < numRowsRational(); i++ )
      {
        assert(_rowTypes[i] == _rangeTypeRational(lhsRational(i), rhsRational(i)));
      }


    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
      _removeRowsReal(perm);

    _invalidateSolution();
}



/// remove all rows with indices in array \p idx of size \p n; an array \p perm of size #numRowsRational() may be
/// passed as buffer memory
template <class R>
void SoPlexBase<R>::removeRowsRational(int idx[], int n, int perm[])
{
    if( perm == 0 )
      {
        DataArray< int > p(numRowsRational());
        _idxToPerm(idx, n, p.get_ptr(), numRowsRational());
        SoPlexBase<R>::removeRowsRational(p.get_ptr());
      }
    else
      {
        _idxToPerm(idx, n, perm, numRowsRational());
        SoPlexBase<R>::removeRowsRational(perm);
      }
}



/// removes rows \p start to \p end including both; an array \p perm of size #numRowsRational() may be passed as
/// buffer memory
template <class R>
void SoPlexBase<R>::removeRowRangeRational(int start, int end, int perm[])
{
    if( perm == 0 )
      {
        DataArray< int > p(numRowsRational());
        _rangeToPerm(start, end, p.get_ptr(), numRowsRational());
        SoPlexBase<R>::removeRowsRational(p.get_ptr());
      }
    else
      {
        _rangeToPerm(start, end, perm, numRowsRational());
        SoPlexBase<R>::removeRowsRational(perm);
      }
}



/// removes column i
template <class R>
void SoPlexBase<R>::removeColRational(int i)
{
    assert(_rationalLP != 0);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->removeCol(i);
    // only swap elements if not the last one was removed
    if( i < _rationalLP->nCols() )
      {
        _colTypes[i] = _colTypes[_rationalLP->nCols()];
        assert(_colTypes[i] == _rangeTypeRational(lowerRational(i), upperRational(i)));
      }
    _colTypes.reSize(_rationalLP->nCols());

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
      _removeColReal(i);

    _invalidateSolution();
}



/// removes all columns with an index \p i such that \p perm[i] < 0; upon completion, \p perm[i] >= 0 indicates the
/// new index where column \p i has been moved to; note that \p perm must point to an array of size at least
/// #numColsRational()
template <class R>
void SoPlexBase<R>::removeColsRational(int perm[])
{
    assert(_rationalLP != 0);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    const int oldsize = numColsRational();
    _rationalLP->removeCols(perm);
    for( int i = 0; i < oldsize; i++ )
      {
        if( perm[i] >= 0 )
          _colTypes[perm[i]] = _colTypes[i];
      }
    _colTypes.reSize(_rationalLP->nCols());
    for( int i = 0; i < numColsRational(); i++ )
      {
        assert(_colTypes[i] == _rangeTypeRational(lowerRational(i), upperRational(i)));
      }

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
      _removeColsReal(perm);

    _invalidateSolution();
}



/// remove all columns with indices in array \p idx of size \p n; an array \p perm of size #numColsRational() may be
/// passed as buffer memory
template <class R>
void SoPlexBase<R>::removeColsRational(int idx[], int n, int perm[])
{
    if( perm == 0 )
      {
        DataArray< int > p(numColsRational());
        _idxToPerm(idx, n, p.get_ptr(), numColsRational());
        SoPlexBase<R>::removeColsRational(p.get_ptr());
      }
    else
      {
        _idxToPerm(idx, n, perm, numColsRational());
        SoPlexBase<R>::removeColsRational(perm);
      }
}



/// removes columns \p start to \p end including both; an array \p perm of size #numColsRational() may be passed as
/// buffer memory
template <class R>
void SoPlexBase<R>::removeColRangeRational(int start, int end, int perm[])
{
    if( perm == 0 )
      {
        DataArray< int > p(numColsRational());
        _rangeToPerm(start, end, p.get_ptr(), numColsRational());
        SoPlexBase<R>::removeColsRational(p.get_ptr());
      }
    else
      {
        _rangeToPerm(start, end, perm, numColsRational());
        SoPlexBase<R>::removeColsRational(perm);
      }
}



/// clears the LP
template <class R>
void SoPlexBase<R>::clearLPRational()
{
    assert(_rationalLP != 0);

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
      return;

    _rationalLP->clear();
    _rationalLUSolver.clear();
    _rowTypes.clear();
    _colTypes.clear();

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_AUTO )
      {
        _realLP->clear();
        _hasBasis = false;
      }

    _invalidateSolution();
}



/// synchronizes rational LP with R LP, i.e., copies R LP to rational LP, if sync mode is manual
template <class R>
void SoPlexBase<R>::syncLPRational()
{
    assert(_isConsistent());

    if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_MANUAL )
      _syncLPRational();
}

/// returns the current solver status
template <class R>
typename SPxSolverBase<R>::Status SoPlexBase<R>::status() const
{
    return _status;
}

/// returns the current basis status
template <class R>
typename SPxBasisBase<R>::SPxStatus SoPlexBase<R>::basisStatus() const
{
    if( !hasBasis() )
      return SPxBasisBase<R>::NO_PROBLEM;
    else if( status() == SPxSolverBase<R>::OPTIMAL || status() == SPxSolverBase<R>::OPTIMAL_UNSCALED_VIOLATIONS )
      return SPxBasisBase<R>::OPTIMAL;
    else if( status() == SPxSolverBase<R>::UNBOUNDED )
      return SPxBasisBase<R>::UNBOUNDED;
    else if( status() == SPxSolverBase<R>::INFEASIBLE )
      return SPxBasisBase<R>::INFEASIBLE;
    else if( hasPrimal() )
      return SPxBasisBase<R>::PRIMAL;
    else if( hasDual() )
      return SPxBasisBase<R>::DUAL;
    else
      return SPxBasisBase<R>::REGULAR;
}


/// returns the objective value if a primal or dual solution is available
template <class R>
R SoPlexBase<R>::objValueReal()
{
    assert(OBJSENSE_MAXIMIZE == 1);
    assert(OBJSENSE_MINIMIZE == -1);

    if( status() == SPxSolverBase<R>::UNBOUNDED )
      return RealParam(SoPlexBase<R>::INFTY) * intParam(SoPlexBase<R>::OBJSENSE);
    else if( status() == SPxSolverBase<R>::INFEASIBLE )
      return -realParam(SoPlexBase<R>::INFTY) * intParam(SoPlexBase<R>::OBJSENSE);
    else if( hasPrimal() || hasDual() )
      {
        _syncRealSolution();
        return _solReal._objVal;
      }
    else
      return 0.0;
}


/// returns the objective value if a primal or dual solution is available
template <class R>
Rational SoPlexBase<R>::objValueRational()
{
  assert(OBJSENSE_MAXIMIZE == 1);
  assert(OBJSENSE_MINIMIZE == -1);

  if( this->status() == SPxSolverBase<R>::UNBOUNDED )
    {
      if( intParam(SoPlexBase<R>::OBJSENSE) == OBJSENSE_MAXIMIZE )
        return _rationalPosInfty;
      else
        return _rationalNegInfty;
  }
  else if( this->status() == SPxSolverBase<R>::INFEASIBLE )
    {
      if( intParam(SoPlexBase<R>::OBJSENSE) == OBJSENSE_MAXIMIZE )
        return _rationalNegInfty;
      else
        return _rationalPosInfty;
    }
  else if( hasPrimal() || hasDual() )
    {
      _syncRationalSolution();
      return _solRational._objVal;
    }
  else
    return _rationalZero;
}

/// Is stored primal solution feasible?
template <class R>
bool SoPlexBase<R>::isPrimalFeasible() const
{
  return (_hasSolReal && _solReal.isPrimalFeasible()) || (_hasSolRational && _solRational.isPrimalFeasible());
}


/// is a primal feasible solution available?
template <class R>
bool SoPlexBase<R>::hasPrimal() const
{
  return _hasSolReal || _hasSolRational;
}

/// is a primal unbounded ray available?
template <class R>
bool SoPlexBase<R>::hasPrimalRay() const
{
  return (_hasSolReal && _solReal.hasPrimalRay()) || (_hasSolRational && _solRational.hasPrimalRay());
}


/// is stored dual solution feasible?
template <class R>
bool SoPlexBase<R>::isDualFeasible() const
{
  return (_hasSolReal && _solReal.isDualFeasible()) || (_hasSolRational && _solRational.isDualFeasible());

}

/// is a dual feasible solution available?
template <class R>
bool SoPlexBase<R>::hasDual() const
{
  return _hasSolReal || _hasSolRational;
}

/// is Farkas proof of infeasibility available?
template <class R>
bool SoPlexBase<R>::hasDualFarkas() const
{
  return (_hasSolReal && _solReal.hasDualFarkas()) || (_hasSolRational && _solRational.hasDualFarkas());
}

/// gets the vector of slack values if available; returns true on success
template <class R>
bool SoPlexBase<R>::getSlacksReal(VectorBase<R>& vector)
{
  if( hasPrimal() && vector.dim() >= numRows() )
    {
      _syncRealSolution();
      _solReal.getSlacks(vector);
      return true;
    }
  else
    return false;
}

/// gets the primal solution vector if available; returns true on success
template <class R>
bool SoPlexBase<R>::getPrimalRational(VectorBase<Rational>& vector)
{
  if( _rationalLP != 0 && hasPrimal() && vector.dim() >= numColsRational() )
    {
      _syncRationalSolution();
      _solRational.getPrimalSol(vector);
      return true;
    }
  else
    return false;
}


/// gets the vector of slack values if available; returns true on success
template <class R>
bool SoPlexBase<R>::getSlacksRational(VectorRational& vector)
{
  if( _rationalLP != 0 && hasPrimal() && vector.dim() >= numRowsRational() )
    {
      _syncRationalSolution();
      _solRational.getSlacks(vector);
      return true;
    }
  else
    return false;
}

template <class R>
bool SoPlexBase<R>::getPrimalRayRational(VectorBase<Rational>& vector)
{
  if( _rationalLP != 0 && hasPrimalRay() && vector.dim() >= numColsRational() )
    {
      _syncRationalSolution();
      _solRational.getPrimalRaySol(vector);
      return true;
    }
  else
    return false;
}



/// gets the dual solution vector if available; returns true on success
template <class R>
bool SoPlexBase<R>::getDualRational(VectorBase<Rational>& vector)
{
  if( _rationalLP != 0 && hasDual() && vector.dim() >= numRowsRational() )
    {
      _syncRationalSolution();
      _solRational.getDualSol(vector);
      return true;
    }
  else
    return false;
}



/// gets the vector of reduced cost values if available; returns true on success
template <class R>
bool SoPlexBase<R>::getRedCostRational(VectorRational& vector)
{
  if( _rationalLP != 0 && hasDual() && vector.dim() >= numColsRational() )
    {
      _syncRationalSolution();
      _solRational.getRedCostSol(vector);
      return true;
    }
  else
    return false;
}

/// gets the Farkas proof if LP is infeasible; returns true on success
template <class R>
bool SoPlexBase<R>::getDualFarkasRational(VectorBase<Rational>& vector)
{
  if( _rationalLP != 0 && hasDualFarkas() && vector.dim() >= numRowsRational() )
    {
      _syncRationalSolution();
      _solRational.getDualFarkasSol(vector);
      return true;
    }
  else
    return false;
}



/// gets violation of bounds; returns true on success
template <class R>
bool SoPlexBase<R>::getBoundViolationRational(Rational& maxviol, Rational& sumviol)
{
  if( !isPrimalFeasible() )
    return false;

  // if we have to synchronize, we do not measure time, because this would affect the solving statistics
  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    _syncLPRational(false);

  _syncRationalSolution();
  VectorRational& primal = _solRational._primal;
  assert(primal.dim() == numColsRational());

  maxviol = 0;
  sumviol = 0;

  for( int i = numColsRational() - 1; i >= 0; i-- )
    {
      Rational viol = lowerRational(i) - primal[i];
      if( viol > 0 )
        {
          sumviol += viol;
          if( viol > maxviol )
            {
              maxviol = viol;
              MSG_DEBUG( std::cout << "increased bound violation for column " << i << ": " << rationalToString(viol)
                         << " lower: " << rationalToString(lowerRational(i))
                         << ", primal: " << rationalToString(primal[i]) << "\n" );
            }
        }

      viol = primal[i] - upperRational(i);
      if( viol > 0 )
        {
          sumviol += viol;
          if( viol > maxviol )
            {
              maxviol = viol;
              MSG_DEBUG( std::cout << "increased bound violation for column " << i << ": " << rationalToString(viol)
                         << " upper: " << rationalToString(upperRational(i))
                         << ", primal: " << rationalToString(primal[i]) << "\n" );
            }
        }
    }

  return true;
}



/// gets violation of constraints; returns true on success
template <class R>
bool SoPlexBase<R>::getRowViolationRational(Rational& maxviol, Rational& sumviol)
{
  if( !isPrimalFeasible() )
    return false;

  // if we have to synchronize, we do not measure time, because this would affect the solving statistics
  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    _syncLPRational(false);

  _syncRationalSolution();
  VectorRational& primal = _solRational._primal;
  assert(primal.dim() == numColsRational());

  DVectorRational activity(numRowsRational());
  _rationalLP->computePrimalActivity(primal, activity);
  maxviol = 0;
  sumviol = 0;

  for( int i = numRowsRational() - 1; i >= 0; i-- )
    {
      Rational viol = lhsRational(i) - activity[i];
      if( viol > 0 )
        {
          sumviol += viol;
          if( viol > maxviol )
            {
              maxviol = viol;
              MSG_DEBUG( std::cout << "increased constraint violation for row " << i << ": " << rationalToString(viol)
                         << " lhs: " << rationalToString(lhsRational(i))
                         << ", activity: " << rationalToString(activity[i]) << "\n" );
            }
        }

      viol = activity[i] - rhsRational(i);
      if( viol > 0 )
        {
          sumviol += viol;
          if( viol > maxviol )
            {
              maxviol = viol;
              MSG_DEBUG( std::cout << "increased constraint violation for row " << i << ": " << rationalToString(viol)
                         << " rhs: " << rationalToString(rhsRational(i))
                         << ", activity: " << rationalToString(activity[i]) << "\n" );
            }
        }
    }

  return true;
}



/// gets violation of reduced costs; returns true on success
template <class R>
bool SoPlexBase<R>::getRedCostViolationRational(Rational& maxviol, Rational& sumviol)
{
  if( !isPrimalFeasible() || !isDualFeasible() )
    return false;

  // if we have to synchronize, we do not measure time, because this would affect the solving statistics
  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    _syncLPRational(false);

  _syncRationalSolution();
  VectorRational& redcost = _solRational._redCost;
  assert(redcost.dim() == numColsRational());

  maxviol = 0;
  sumviol = 0;

  for( int c = numCols() - 1; c >= 0; c-- )
    {
      assert(!_hasBasis || basisColStatus(c) != SPxSolverBase<R>::UNDEFINED);

      if( _colTypes[c] == RANGETYPE_FIXED )
        {
          assert(lowerRational(c) == upperRational(c));
          continue;
        }

      assert(!_hasBasis || basisColStatus(c) != SPxSolverBase<R>::ON_LOWER || _solRational._primal[c] == lowerRational(c));
      assert(!_hasBasis || basisColStatus(c) != SPxSolverBase<R>::ON_UPPER || _solRational._primal[c] == upperRational(c));
      assert(!_hasBasis || basisColStatus(c) != SPxSolverBase<R>::FIXED || _solRational._primal[c] == lowerRational(c));
      assert(!_hasBasis || basisColStatus(c) != SPxSolverBase<R>::FIXED || _solRational._primal[c] == upperRational(c));

      if( intParam(SoPlexBase<R>::OBJSENSE) == OBJSENSE_MINIMIZE )
        {
          if( _solRational._primal[c] != upperRational(c) && redcost[c] < 0 )
            {
              sumviol += -redcost[c];
              if( redcost[c] < -maxviol )
                {
                  MSG_DEBUG( std::cout << "increased reduced cost violation for column " << c << " not on upper bound: " << rationalToString(-redcost[c]) << "\n" );
                  maxviol = -redcost[c];
                }
            }
          if( _solRational._primal[c] != lowerRational(c) && redcost[c] > 0 )
            {
              sumviol += redcost[c];
              if( redcost[c] > maxviol )
                {
                  MSG_DEBUG( std::cout << "increased reduced cost violation for column " << c << " not on lower bound: " << rationalToString(redcost[c]) << "\n" );
                  maxviol = redcost[c];
                }
            }
        }
      else
        {
          if( _solRational._primal[c] != upperRational(c) && redcost[c] > 0 )
            {
              sumviol += redcost[c];
              if( redcost[c] > maxviol )
                {
                  MSG_DEBUG( std::cout << "increased reduced cost violation for column " << c << " not on upper bound: " << rationalToString(redcost[c]) << "\n" );
                  maxviol = redcost[c];
                }
            }
          if( _solRational._primal[c] != lowerRational(c) && redcost[c] < 0 )
            {
              sumviol += -redcost[c];
              if( redcost[c] < -maxviol )
                {
                  MSG_DEBUG( std::cout << "increased reduced cost violation for column " << c << " not on lower bound: " << rationalToString(-redcost[c]) << "\n" );
                  maxviol = -redcost[c];
                }
            }
        }
    }

  return true;
}



/// gets violation of dual multipliers; returns true on success
template <class R>
bool SoPlexBase<R>::getDualViolationRational(Rational& maxviol, Rational& sumviol)
{
  if( !isDualFeasible() || !isPrimalFeasible() )
    return false;

  // if we have to synchronize, we do not measure time, because this would affect the solving statistics
  if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL )
    _syncLPRational(false);

  _syncRationalSolution();
  VectorRational& dual = _solRational._dual;
  assert(dual.dim() == numRowsRational());

  maxviol = 0;
  sumviol = 0;

  for( int r = numRows() - 1; r >= 0; r-- )
    {
      assert(!_hasBasis || basisRowStatus(r) != SPxSolverBase<R>::UNDEFINED);

      if( _rowTypes[r] == RANGETYPE_FIXED )
        {
          assert(lhsRational(r) == rhsRational(r));
          continue;
        }

      assert(!_hasBasis || basisRowStatus(r) != SPxSolverBase<R>::ON_LOWER || _solRational._slacks[r] <= lhsRational(r) + _rationalFeastol);
      assert(!_hasBasis || basisRowStatus(r) != SPxSolverBase<R>::ON_UPPER || _solRational._slacks[r] >= rhsRational(r) - _rationalFeastol);
      assert(!_hasBasis || basisRowStatus(r) != SPxSolverBase<R>::FIXED || _solRational._slacks[r] <= lhsRational(r) + _rationalFeastol);
      assert(!_hasBasis || basisRowStatus(r) != SPxSolverBase<R>::FIXED || _solRational._slacks[r] >= rhsRational(r) - _rationalFeastol);

      if( intParam(SoPlexBase<R>::OBJSENSE) == OBJSENSE_MINIMIZE )
        {
          if( _solRational._slacks[r] < rhsRational(r) - _rationalFeastol && dual[r] < 0 )
            {
              sumviol += -dual[r];
              if( dual[r] < -maxviol )
                {
                  MSG_DEBUG( std::cout << "increased dual violation for row " << r << " not on upper bound: " << rationalToString(-dual[r])
                             << " (slack = " << rationalToString(_solRational._slacks[r])
                             << ", status = " << basisRowStatus(r)
                             << ", lhs = " << rationalToString(lhsRational(r))
                             << ", rhs = " << rationalToString(rhsRational(r)) << ")\n" );
                  maxviol = -dual[r];
                }
            }
          if( _solRational._slacks[r] > lhsRational(r) + _rationalFeastol && dual[r] > 0 )
            {
              sumviol += dual[r];
              if( dual[r] > maxviol )
                {
                  MSG_DEBUG( std::cout << "increased dual violation for row " << r << " not on lower bound: " << rationalToString(dual[r])
                             << " (slack = " << rationalToString(_solRational._slacks[r])
                             << ", status = " << basisRowStatus(r)
                             << ", lhs = " << rationalToString(lhsRational(r))
                             << ", rhs = " << rationalToString(rhsRational(r)) << ")\n" );
                  maxviol = dual[r];
                }
            }
        }
      else
        {
          if( _solRational._slacks[r] < rhsRational(r) - _rationalFeastol && dual[r] > 0 )
            {
              sumviol += dual[r];
              if( dual[r] > maxviol )
                {
                  MSG_DEBUG( std::cout << "increased dual violation for row " << r << " not on upper bound: " << rationalToString(dual[r])
                             << " (slack = " << rationalToString(_solRational._slacks[r])
                             << ", status = " << basisRowStatus(r)
                             << ", lhs = " << rationalToString(lhsRational(r))
                             << ", rhs = " << rationalToString(rhsRational(r)) << ")\n" );
                  maxviol = dual[r];
                }
            }
          if( _solRational._slacks[r] > lhsRational(r) + _rationalFeastol && dual[r] < 0 )
            {
              sumviol += -dual[r];
              if( dual[r] < -maxviol )
                {
                  MSG_DEBUG( std::cout << "increased dual violation for row " << r << " not on lower bound: " << rationalToString(-dual[r])
                             << " (slack = " << rationalToString(_solRational._slacks[r])
                             << ", status = " << basisRowStatus(r)
                             << ", lhs = " << rationalToString(lhsRational(r))
                             << ", rhs = " << rationalToString(rhsRational(r)) << ")\n" );
                  maxviol = -dual[r];
                }
            }
        }
    }

  return true;
}


  /// get size of primal solution
  template <class R>
  int SoPlexBase<R>::totalSizePrimalRational(const int base)
  {
    if( hasPrimal() || hasPrimalRay() )
      {
        _syncRationalSolution();
        return _solRational.totalSizePrimal(base);
      }
    else
      return 0;
  }



  /// get size of dual solution
  template <class R>
  int SoPlexBase<R>::totalSizeDualRational(const int base)
  {
    if( hasDual() || hasDualFarkas() )
      {
        _syncRationalSolution();
        return _solRational.totalSizeDual(base);
      }
    else
      return 0;
  }



  /// get size of least common multiple of denominators in primal solution
  template <class R>
  int SoPlexBase<R>::dlcmSizePrimalRational(const int base)
  {
    if( hasPrimal() || hasPrimalRay() )
      {
        _syncRationalSolution();
        return _solRational.dlcmSizePrimal(base);
      }
    else
      return 0;
  }



  /// get size of least common multiple of denominators in dual solution
  template <class R>
  int SoPlexBase<R>::dlcmSizeDualRational(const int base)
  {
    if( hasDual() || hasDualFarkas() )
      {
        _syncRationalSolution();
        return _solRational.dlcmSizeDual(base);
      }
    else
      return 0;
  }



  /// get size of largest denominator in primal solution
  template <class R>
  int SoPlexBase<R>::dmaxSizePrimalRational(const int base)
  {
    if( hasPrimal() || hasPrimalRay() )
      {
        _syncRationalSolution();
        return _solRational.dmaxSizePrimal(base);
      }
    else
      return 0;
  }



  /// get size of largest denominator in dual solution
  template <class R>
	int SoPlexBase<R>::dmaxSizeDualRational(const int base)
  {
    if( hasDual() || hasDualFarkas() )
      {
        _syncRationalSolution();
        return _solRational.dmaxSizeDual(base);
      }
    else
      return 0;
  }

  /// is an advanced starting basis available?
  template <class R>
	bool SoPlexBase<R>::hasBasis() const
  {
    return _hasBasis;
  }



  /// returns basis status for a single row
  template <class R>
  typename SPxSolverBase<R>::VarStatus SoPlexBase<R>::basisRowStatus(int row) const
  {
    assert(row >= 0);
    assert(row < numRows());

    // if no basis is available, return slack basis; if index is out of range, return basic status as for a newly
    // added row
    if( !hasBasis() || row < 0 || row >= numRows() )
      return SPxSolverBase<R>::BASIC;
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
  template <class R>
  typename SPxSolverBase<R>::VarStatus SoPlexBase<R>::basisColStatus(int col) const
  {
    assert(col >= 0);
    assert(col < numCols());

    // if index is out of range, return nonbasic status as for a newly added unbounded column
    if( col < 0 || col >= numCols() )
      {
        return SPxSolverBase<R>::ZERO;
      }
    // if no basis is available, return slack basis
    else if( !hasBasis() )
      {
        if( lowerReal(col) > -realParam(SoPlexBase<R>::INFTY) )
          return SPxSolverBase<R>::ON_LOWER;
        else if( upperReal(col) < realParam(SoPlexBase<R>::INFTY) )
          return SPxSolverBase<R>::ON_UPPER;
        else
          return SPxSolverBase<R>::ZERO;
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
  template <class R>
  void SoPlexBase<R>::getBasis(typename SPxSolverBase<R>::VarStatus rows[], typename SPxSolverBase<R>::VarStatus cols[]) const
  {
    // if no basis is available, return slack basis
    if( !hasBasis() )
      {
        for( int i = numRows() - 1; i >= 0; i-- )
          rows[i] = SPxSolverBase<R>::BASIC;

        for( int i = numCols() - 1; i >= 0; i-- )
          {
            if( lowerReal(i) > -realParam(SoPlexBase<R>::INFTY) )
              cols[i] = SPxSolverBase<R>::ON_LOWER;
            else if( upperReal(i) < realParam(SoPlexBase<R>::INFTY) )
              cols[i] = SPxSolverBase<R>::ON_UPPER;
            else
              cols[i] = SPxSolverBase<R>::ZERO;
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
        assert(numRows() == _basisStatusRows.size());
        assert(numCols() == _basisStatusCols.size());

        for( int i = numRows() - 1; i >= 0; i-- )
          rows[i] = _basisStatusRows[i];

        for( int i = numCols() - 1; i >= 0; i-- )
          cols[i] = _basisStatusCols[i];
      }
  }



  /// returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m
  template <class R>
  void SoPlexBase<R>::getBasisInd(int* bind) const
  {
    // if no basis is available, return slack basis
    if( !hasBasis() )
      {
        for( int i = 0; i < numRows(); ++i )
          bind[i] = -1 - i;
      }
    // if the real LP is not loaded, the basis is stored in the basis arrays of this class
    else if( !_isRealLPLoaded )
      {
        int k = 0;

        assert(numRows() == _basisStatusRows.size());
        assert(numCols() == _basisStatusCols.size());

        for( int i = 0; i < numRows(); ++i )
          {
            if( _basisStatusRows[i] == SPxSolverBase<R>::BASIC )
              {
                bind[k] = -1 - i;
                k++;
              }
          }

        for( int j = 0; j < numCols(); ++j )
          {
            if( _basisStatusCols[j] == SPxSolverBase<R>::BASIC )
              {
                bind[k] = j;
                k++;
              }
          }

        assert(k == numRows());
      }
    // if the real LP is loaded, the basis is stored in the solver and we need to distinguish between column and row
    // representation; ask the solver itself which representation it has, since the REPRESENTATION parameter of this
    // class might be set to automatic
    else if( _solver.rep() == SPxSolverBase<R>::COLUMN )
      {
        for( int i = 0; i < numRows(); ++i )
          {
            SPxId id = _solver.basis().baseId(i);
            bind[i] = (id.isSPxColId() ? _solver.number(id) : - 1 - _solver.number(id));
          }
      }
    // for row representation, return the complement of the row basis; for this, we need to loop through all rows and columns
    else
      {
        assert(_solver.rep() == SPxSolverBase<R>::ROW);

        int k = 0;

        for( int i = 0; i < numRows(); ++i )
          {
            if( !_solver.isRowBasic(i) )
              {
                bind[k] = -1 - i;
                k++;
              }
          }

        for( int j = 0; j < numCols(); ++j )
          {
            if( !_solver.isColBasic(j) )
            {
               bind[k] = j;
               k++;
            }
         }

         assert(k == numRows());
      }
   }


   /** compute one of several matrix metrics based on the diagonal of the LU factorization
    *  type = 0: max/min ratio
    *  type = 1: trace of U (sum of diagonal elements)
    *  type = 2: determinant (product of diagonal elements)
    */
  template <class R>
  bool SoPlexBase<R>::getBasisMetric(R& condition, int type)
   {
      _ensureRealLPLoaded();
      if( !_isRealLPLoaded )
         return false;

    if( _solver.basis().status() == SPxBasisBase<R>::NO_PROBLEM )
      {
      return false;
      }

      condition = _solver.getBasisMetric(type);

    return true;
  }

  /// computes an estimated condition number for the current basis matrix using the power method; returns true on success
  template <class R>
	bool SoPlexBase<R>::getEstimatedCondition(R& condition)
  {
    _ensureRealLPLoaded();
    if( !_isRealLPLoaded )
      return false;

    if( _solver.basis().status() == SPxBasisBase<R>::NO_PROBLEM )
      return false;

    condition = _solver.basis().getEstimatedCondition();

    return true;
  }

  /// computes the exact condition number for the current basis matrix using the power method; returns true on success
  template <class R>
	bool SoPlexBase<R>::getExactCondition(R& condition)
  {
    _ensureRealLPLoaded();
    if( !_isRealLPLoaded )
      return false;

    if( _solver.basis().status() == SPxBasisBase<R>::NO_PROBLEM )
      return false;

    condition = _solver.basis().getExactCondition();

    return true;
  }

  /// computes row r of basis inverse; returns true on success
  template <class R>
	bool SoPlexBase<R>::getBasisInverseRowReal(int r, R* coef, int* inds, int* ninds, bool unscale)
  {
    assert(r >= 0);
    assert(r < numRows());
    assert(coef != 0);

    if( !hasBasis() || r < 0 || r >= numRows() )
      return false;

    _ensureRealLPLoaded();

    if( !_isRealLPLoaded )
      return false;

    // we need to distinguish between column and row representation; ask the solver itself which representation it
    // has, since the REPRESENTATION parameter of this class might be set to automatic
    if( _solver.rep() == SPxSolverBase<R>::COLUMN )
      {
        int idx;
        SSVectorBase<R> x(numRows());
        try
          {
            /* unscaling required? */
            if( unscale && _solver.isScaled())
              {
                /* for information on the unscaling procedure see spxscaler.h */

                int scaleExp;
                DSVector rhs(_solver.unitVector(r));

                if( _solver.basis().baseId(r).isSPxColId() )
                  scaleExp = _scaler->getColScaleExp(_solver.number(_solver.basis().baseId(r)));
                else
                  scaleExp = - _scaler->getRowScaleExp(_solver.number(_solver.basis().baseId(r)));

                rhs *= spxLdexp(1.0, scaleExp);

                _solver.basis().coSolve(x, rhs);
                x.setup();
                int size = x.size();

                for( int i = 0; i < size; i++ )
                  {
                    scaleExp = _scaler->getRowScaleExp(x.index(i));
                    x.scaleValue(x.index(i), scaleExp);
                  }
              }
            else
              {
                _solver.basis().coSolve(x, _solver.unitVector(r));
              }
          }
        catch( const SPxException& E )
          {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
          }
        // copy sparse data to dense result vector based on coef array
        if( ninds != 0 && inds != 0 )
          {
            // during solving SoPlexBase may have destroyed the sparsity structure so we need to restore it
            x.setup();
            *ninds = x.size();
            for( int i = 0; i < *ninds; ++i )
              {
                idx = x.index(i);
                coef[idx] = x[idx];
                // set sparsity pattern of coef array
                inds[i] = idx;
              }
          }
        else
          {
            VectorBase<R> y(numRows(), coef);
            y = x;
            if( ninds != NULL )
              *ninds = -1;
          }
      }
    else
      {
        assert(_solver.rep() == SPxSolverBase<R>::ROW);

        // @todo should rhs be a reference?
        DSVector rhs(numCols());
        SSVector y(numCols());
        int* bind = 0;
        int index;

        // get ordering of column basis matrix
        spx_alloc(bind, numRows());
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
            assert(index < numRows());
            assert(!_solver.isRowBasic(index));

            // get row vector
            rhs = _solver.rowVector(index);
            rhs *= -1.0;

            if( unscale && _solver.isScaled() )
              {
                for( int i = 0; i < rhs.size(); ++i)
                  rhs.value(i) = spxLdexp(rhs.value(i), -_scaler->getRowScaleExp(index));
              }
          }
        // r corresponds to a column vector
        else
          {
            // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index < numCols());
            assert(!_solver.isColBasic(index));

            // get unit vector
            rhs = UnitVectorReal(index);

            if( unscale && _solver.isScaled() )
              rhs *= spxLdexp(1.0, _scaler->getColScaleExp(index));
          }

        // solve system "y B = rhs", where B is the row basis matrix
        try
          {
            _solver.basis().solve(y, rhs);
          }
        catch( const SPxException& E )
          {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
          }

        // initialize result vector x as zero
        memset(coef, 0, (unsigned int)numRows() * sizeof(Real));

        // add nonzero entries
        for( int i = 0; i < numCols(); ++i )
          {
            SPxId id = _solver.basis().baseId(i);

            if( id.isSPxRowId() )
              {
                assert(_solver.number(id) >= 0);
                assert(_solver.number(id) < numRows());
                assert(bind[r] >= 0 || _solver.number(id) != index);

                int rowindex = _solver.number(id);
                coef[rowindex] = y[i];

                if( unscale && _solver.isScaled() )
                  coef[rowindex] = spxLdexp(y[i], _scaler->getRowScaleExp(rowindex));
              }
          }

        // if r corresponds to a row vector, we have to add a 1 at position r
        if( bind[r] < 0 )
          {
            assert(coef[index] == 0.0);
            coef[index] = 1.0;
          }

        // @todo implement returning of sparsity information like in column wise case
        if( ninds != NULL)
          *ninds = -1;

        // free memory
        spx_free(bind);
      }

    return true;
  }



  /// computes column c of basis inverse; returns true on success
  /// @todo does not work correctly for the row representation
  template <class R>
	bool SoPlexBase<R>::getBasisInverseColReal(int c, R* coef, int* inds, int* ninds, bool unscale)
  {
    assert(c >= 0);
    assert(c < numRows());
    assert(coef != 0);

    if( !hasBasis() || c < 0 || c >= numRows() )
      return false;

    _ensureRealLPLoaded();

    if( !_isRealLPLoaded )
      return false;

    // we need to distinguish between column and row representation; ask the solver itself which representation it
    // has, since the REPRESENTATION parameter of this class might be set to automatic
    if( _solver.rep() == SPxSolverBase<R>::COLUMN )
      {
        int idx;
        SSVectorBase<R> x(numRows());
        try
          {
            /* unscaling required? */
            if( unscale && _solver.isScaled())
              {
                /* for information on the unscaling procedure see spxscaler.h */

                int scaleExp =_scaler->getRowScaleExp(c);
                DSVector rhs(_solver.unitVector(c));
                rhs *= spxLdexp(1.0, scaleExp);

                _solver.basis().solve(x, rhs);

                x.setup();
                int size = x.size();

                for( int i = 0; i < size; i++ )
                  {
                    if( _solver.basis().baseId(x.index(i)).isSPxColId() )
                      {
                        idx = _solver.number(_solver.basis().baseId(x.index(i)));
                        scaleExp = _scaler->getColScaleExp(idx);
                        x.scaleValue(x.index(i), scaleExp);
                      }
                    else
                      {
                        idx = _solver.number(_solver.basis().baseId(x.index(i)));
                        scaleExp = - _scaler->getRowScaleExp(idx);
                        x.scaleValue(x.index(i), scaleExp);
                      }
                  }
              }
            else
              {
                _solver.basis().solve(x, _solver.unitVector(c));
              }
          }
        catch( const SPxException& E )
          {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
          }
        // copy sparse data to dense result vector based on coef array
        if( ninds != 0 && inds != 0 )
          {
            // SoPlexBase may have destroyed the sparsity structure so we need to restore it
            x.setup();
            *ninds = x.size();
            for( int i = 0; i < *ninds; ++i )
              {
                idx = x.index(i);
                coef[idx] = x[idx];
                // set sparsity pattern of coef array
                inds[i] = idx;
              }
          }
        else
          {
            VectorBase<R> y(numRows(), coef);
            y = x;
            if( ninds != 0 )
              *ninds = -1;
          }
      }
    else
      {
        assert(_solver.rep() == SPxSolverBase<R>::ROW);

        // @todo should rhs be a reference?
        DSVectorBase<R> rhs(numCols());
        SSVectorBase<R> y(numCols());
        int* bind = 0;
        int index;

        // get ordering of column basis matrix
        spx_alloc(bind, numRows());
        getBasisInd(bind);

        // get vector corresponding to requested index c
        index = bind[c];

        // c corresponds to a row vector
        if( index < 0 )
          {
            // transform index to actual row index
            index = -index - 1;

            // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index >= 0);
            assert(index < numRows());
            assert(!_solver.isRowBasic(index));

            // get row vector
            rhs = _solver.rowVector(index);
            rhs *= -1.0;
          }
        // c corresponds to a column vector
        else
          {
            // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
            assert(index < numCols());
            assert(!_solver.isColBasic(index));

            // get unit vector
            rhs = UnitVectorReal(index);
          }

        // solve system "y B = rhs", where B is the row basis matrix
        try
          {
            /* unscaling required? */
            if( unscale && _solver.isScaled() )
              {
                int size = rhs.size();
                int scaleExp;

                for( int i = 0; i < size; i++ )
                  {
                    scaleExp = _scaler->getColScaleExp(i);
                    rhs.value(i) *= spxLdexp(1.0, scaleExp);
                  }

                _solver.basis().coSolve(y, rhs);

                int rowIdx;
                size = y.size();

                for( int i = 0; i < size; i++ )
                  {
                    assert(_solver.basis().baseId(y.index(i)).isSPxRowId());
                    rowIdx = _solver.basis().baseId(y.index(i)).getIdx();
                    scaleExp = _scaler->getRowScaleExp(rowIdx);
                    y.setValue(i, y.value(i) * spxLdexp(1.0, scaleExp));
                  }
              }
            else
              {
                _solver.basis().coSolve(y, rhs);
              }
          }
        catch( const SPxException& E )
          {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while computing basis inverse row.\n" );
            return false;
          }

        // initialize result vector x as zero
        memset(coef, 0, (unsigned int)numRows() * sizeof(Real));

        // add nonzero entries
        for( int i = 0; i < numCols(); ++i )
          {
            SPxId id = _solver.basis().baseId(i);

            if( id.isSPxRowId() )
              {
                assert(_solver.number(id) >= 0);
                assert(_solver.number(id) < numRows());
                assert(bind[c] >= 0 || _solver.number(id) != index);

                coef[_solver.number(id)] = y[i];
              }
          }

        // if c corresponds to a row vector, we have to add a 1 at position c
        if( bind[c] < 0 )
          {
            assert(coef[index] == 0.0);
            coef[index] = 1.0;
          }

        // @todo implement returning of sparsity information like in column wise case
        if( ninds != NULL)
          *ninds = -1;

        // free memory
        spx_free(bind);
      }

    return true;
  }



  /// computes dense solution of basis matrix B * sol = rhs; returns true on success
  template <class R>
	bool SoPlexBase<R>::getBasisInverseTimesVecReal(Real* rhs, R* sol, bool unscale)
  {
    VectorBase<R> v(numRows(), rhs);
    VectorBase<R> x(numRows(), sol);

    if( !hasBasis() )
      return false;

    _ensureRealLPLoaded();

    if( !_isRealLPLoaded )
      return false;

    // we need to distinguish between column and row representation; ask the solver itself which representation it
    // has, since the REPRESENTATION parameter of this class might be set to automatic; in the column case we can use
    // the existing factorization
    if( _solver.rep() == SPxSolverBase<R>::COLUMN )
      {
        // solve system "x = B^-1 * v"
        try
          {
            /* unscaling required? */
            if( unscale && _solver.isScaled())
              {
                /* for information on the unscaling procedure see spxscaler.h */
                int scaleExp;
                int idx;

                for( int i = 0; i < v.dim(); ++i)
                  {
                    if( isNotZero(v[i]) )
                      {
                        scaleExp =_scaler->getRowScaleExp(i);
                        v[i] = spxLdexp(v[i], scaleExp);
                      }
                  }

                _solver.basis().solve(x, v);

                for( int i = 0; i < x.dim(); i++ )
                  {
                    if( isNotZero(x[i]) )
                      {
                        idx = _solver.number(_solver.basis().baseId(i));
                        if( _solver.basis().baseId(i).isSPxColId() )
                          scaleExp = _scaler->getColScaleExp(idx);
                        else
                          scaleExp = - _scaler->getRowScaleExp(idx);
                        x[i] = spxLdexp(x[i], scaleExp);
                      }
                  }
              }
            else
              {
                _solver.basis().solve(x, v);
              }
          }
        catch( const SPxException& E )
          {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while solving with basis matrix.\n" );
            return false;
          }
      }
    else
      {
        assert(_solver.rep() == SPxSolverBase<R>::ROW);

        DSVectorBase<R> rowrhs(numCols());
        SSVectorBase<R> y(numCols());
        int* bind = 0;

        bool adaptScaling = unscale && _realLP->isScaled();
        int scaleExp;
        int idx;

        // get ordering of column basis matrix
        spx_alloc(bind, numRows());
        getBasisInd(bind);

        // fill right-hand side for row-based system
        for( int i = 0; i < numCols(); ++i )
          {
            SPxId id = _solver.basis().baseId(i);

            if( id.isSPxRowId() )
              {
                assert(_solver.number(id) >= 0);
                assert(_solver.number(id) < numRows());

                if( adaptScaling )
                  {
                    idx = _solver.number(id);
                    scaleExp = _scaler->getRowScaleExp(idx);
                    rowrhs.add(i, spxLdexp(v[idx], scaleExp));
                  }
                else
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
        catch( const SPxException& E )
          {
            MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while solving with basis matrix.\n" );
            return false;
          }

        // fill result w.r.t. order given by bind
        for( int i = 0; i < numRows(); ++i )
          {
            int index;

            index = bind[i];

            if( index < 0 )
              {
                index = -index-1;

                // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
                assert(index >= 0);
                assert(index < numRows());
                assert(!_solver.isRowBasic(index));

                x[i] = v[index] - (rowVectorRealInternal(index) * Vector(numCols(), y.get_ptr()));

                if( adaptScaling )
                  {
                    scaleExp = -_scaler->getRowScaleExp(index);
                    x[i] = spxLdexp(x[i], scaleExp);
                  }
              }
            else
              {
                // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
                assert(index >= 0);
                assert(index < numCols());
                assert(!_solver.isColBasic(index));

                if( adaptScaling )
                  {
                    scaleExp = _scaler->getColScaleExp(index);
                    x[i] = spxLdexp(y[index], scaleExp);
                  }
                else
                  x[i] = y[index];
              }
          }

        // free memory
        spx_free(bind);
      }
    return true;
  }



  /// multiply with basis matrix; B * vec (inplace)
  template <class R>
	bool SoPlexBase<R>::multBasis(Real* vec, bool unscale)
  {
    if( !hasBasis() )
      return false;

    _ensureRealLPLoaded();

    if( !_isRealLPLoaded )
      return false;

    if( _solver.rep() == SPxSolverBase<R>::COLUMN )
      {
        int basisdim = numRows();

        // create Vector from input values
        Vector x(basisdim, vec);

        if( unscale && _solver.isScaled() )
          {
            /* for information on the unscaling procedure see spxscaler.h */

            int scaleExp;
            for( int i = 0; i < basisdim; ++i)
              {
                if( isNotZero(vec[i]) )
                  {
                    if( _solver.basis().baseId(i).isSPxColId() )
                      scaleExp = - _scaler->getColScaleExp(_solver.number(_solver.basis().baseId(i)));
                    else
                      scaleExp = _scaler->getRowScaleExp(_solver.number(_solver.basis().baseId(i)));

                    vec[i] = spxLdexp(vec[i], scaleExp);
                  }
              }
            _solver.basis().multBaseWith(x);
            for( int i = 0; i < basisdim; ++i)
              {
                scaleExp = _scaler->getRowScaleExp(i);
                vec[i] = spxLdexp(vec[i], -scaleExp);
              }
          }
        else
          _solver.basis().multBaseWith(x);
      }
    else
      {
        int colbasisdim = numRows();

        DSVector y(colbasisdim);

        y.clear();

        // create Vector from input values
        Vector x(colbasisdim, vec);

        int* bind = 0;
        int index;

        // get ordering of column basis matrix
        spx_alloc(bind, colbasisdim);
        getBasisInd(bind);

        // temporarily create the column basis and multiply every column with x
        for( int i = 0; i < colbasisdim; ++i)
          {
            if( isNotZero(x[i]) )
              {
                // get vector corresponding to requested index i
                index = bind[i];
                // r corresponds to a row vector
                if( index < 0 )
                  {
                    // transform index to actual row index
                    index = -index - 1;

                    // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
                    assert(index >= 0);
                    assert(index < numRows());
                    assert(!_solver.isRowBasic(index));

                    y.add(x[i] * UnitVectorReal(index));
                  }
                // r corresponds to a column vector
                else
                  {
                    // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
                    assert(index < numCols());
                    assert(!_solver.isColBasic(index));

                    if( unscale && _solver.isScaled() )
                      {
                        DSVector col;
                        _solver.getColVectorUnscaled(index, col);
                        y.add(x[i] * col);
                      }

                    y.add(x[i] * _solver.colVector(index));
                  }
              }
          }
        spx_free(bind);
        x = y;
      }

    return true;
  }



  /// multiply with transpose of basis matrix; vec * B^T (inplace)
  template <class R>
	bool SoPlexBase<R>::multBasisTranspose(Real* vec, bool unscale)
  {
    if( !hasBasis() )
      return false;

    _ensureRealLPLoaded();

    if( !_isRealLPLoaded )
      return false;

    if( _solver.rep() == SPxSolverBase<R>::COLUMN )
      {
        int basisdim = numRows();

        // create Vector from input values
        Vector x(basisdim, vec);

        if( unscale && _solver.isScaled() )
          {
            /* for information on the unscaling procedure see spxscaler.h */

            int scaleExp;
            for( int i = 0; i < basisdim; ++i)
              {
                if( isNotZero(vec[i]) )
                  {
                    scaleExp = - _scaler->getRowScaleExp(i);
                    vec[i] = spxLdexp(vec[i], scaleExp);
                  }
              }
            _solver.basis().multWithBase(x);
            for( int i = 0; i < basisdim; ++i)
              {
                if( isNotZero(vec[i]) )
                  {
                    if( _solver.basis().baseId(i).isSPxColId() )
                      scaleExp = - _scaler->getColScaleExp(_solver.number(_solver.basis().baseId(i)));
                    else
                      scaleExp = _scaler->getRowScaleExp(_solver.number(_solver.basis().baseId(i)));

                    vec[i] = spxLdexp(vec[i], scaleExp);
                  }
              }
          }
        else
          _solver.basis().multWithBase(x);
      }
    else
      {
        int colbasisdim = numRows();

        DSVector y(colbasisdim);

        // create Vector from input values
        Vector x(colbasisdim, vec);

        int* bind = 0;
        int index;

        // get ordering of column basis matrix
        spx_alloc(bind, colbasisdim);
        getBasisInd(bind);

        // temporarily create the column basis and multiply every column with x
        for( int i = 0; i < colbasisdim; ++i)
          {
            // get vector corresponding to requested index i
            index = bind[i];
            // r corresponds to a row vector
            if( index < 0 )
              {
                // transform index to actual row index
                index = -index - 1;

                // should be a valid row index and in the column basis matrix, i.e., not basic w.r.t. row representation
                assert(index >= 0);
                assert(index < numRows());
                assert(!_solver.isRowBasic(index));

                y.add(i, x * UnitVectorReal(index));
              }
            // r corresponds to a column vector
            else
              {
                // should be a valid column index and in the column basis matrix, i.e., not basic w.r.t. row representation
                assert(index < numCols());
                assert(!_solver.isColBasic(index));

                if( unscale && _solver.isScaled() )
                  {
                    DSVector col;
                    _solver.getColVectorUnscaled(index, col);
                    y.add(i, x * col);
                  }
                else
                  y.add(i, x * _solver.colVector(index));
              }
          }
        spx_free(bind);
        x = y;
      }

    return true;
  }



  /// compute rational basis inverse; returns true on success
  template <class R>
	bool SoPlexBase<R>::computeBasisInverseRational()
  {
    if( !hasBasis() )
      {
        _rationalLUSolver.clear();
        assert(_rationalLUSolver.status() == SLinSolverRational::UNLOADED);
        return false;
      }

    if( _rationalLUSolver.status() == SLinSolverRational::UNLOADED
        || _rationalLUSolver.status() == SLinSolverRational::TIME )
      {
        _rationalLUSolverBind.reSize(numRowsRational());
        getBasisInd(_rationalLUSolverBind.get_ptr());
        _computeBasisInverseRational();
      }

    if( _rationalLUSolver.status() == SLinSolverRational::OK )
      return true;

    return false;
  }



  /// gets an array of indices for the columns of the rational basis matrix; bind[i] >= 0 means that the i-th column of
  /// the basis matrix contains variable bind[i]; bind[i] < 0 means that the i-th column of the basis matrix contains
  /// the slack variable for row -bind[i]-1; performs rational factorization if not available; returns true on success
  template <class R>
	bool SoPlexBase<R>::getBasisIndRational(DataArray<int>& bind)
  {
    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      computeBasisInverseRational();

    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      return false;

    bind = _rationalLUSolverBind;
    assert(bind.size() == numRowsRational());
    return true;
  }



  /// computes row r of basis inverse; performs rational factorization if not available; returns true on success
  template <class R>
	bool SoPlexBase<R>::getBasisInverseRowRational(const int r, SSVectorRational& vec)
  {
    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      computeBasisInverseRational();

    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      return false;

    try
      {
        vec.reDim(numRowsRational());
        _rationalLUSolver.solveLeft(vec, *_unitVectorRational(r));
      }
    catch( const SPxException& E )
      {
        MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while computing rational basis inverse row.\n" );
        return false;
      }
    return true;
  }

  /// computes column c of basis inverse; performs rational factorization if not available; returns true on success
  template <class R>
	bool SoPlexBase<R>::getBasisInverseColRational(const int c, SSVectorRational& vec)
  {
    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      computeBasisInverseRational();

    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      return false;

    try
      {
        vec.reDim(numRowsRational());
        _rationalLUSolver.solveRight(vec, *_unitVectorRational(c));
      }
    catch( const SPxException& E )
      {
        MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> while computing rational basis inverse column.\n" );
        return false;
      }
    return true;
  }



  /// computes solution of basis matrix B * sol = rhs; performs rational factorization if not available; returns true
  /// on success
  template <class R>
	bool SoPlexBase<R>::getBasisInverseTimesVecRational(const SVectorRational& rhs, SSVectorRational& sol)
  {
    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      computeBasisInverseRational();

    if( _rationalLUSolver.status() != SLinSolverRational::OK )
      return false;

    try
      {
        sol.reDim(numRowsRational());
        _rationalLUSolver.solveRight(sol, rhs);
      }
    catch( const SPxException& E )
      {
        MSG_INFO1( spxout, spxout << "Caught exception <" << E.what() << "> during right solve with rational basis inverse.\n" );
        return false;
      }
    return true;
  }



  /// sets starting basis via arrays of statuses
  template <class R>
  void SoPlexBase<R>::setBasis(const typename SPxSolverBase<R>::VarStatus rows[], const typename SPxSolverBase<R>::VarStatus cols[])
  {
    _rationalLUSolver.clear();

    if( _isRealLPLoaded )
      {
        assert(numRows() == _solver.nRows());
        assert(numCols() == _solver.nCols());

        _solver.setBasis(rows, cols);
        _hasBasis = (_solver.basis().status() > SPxBasisBase<R>::NO_PROBLEM);
      }
    else
      {
        _basisStatusRows.reSize(numRows());
        _basisStatusCols.reSize(numCols());

        for( int i = numRows() - 1; i >= 0; i-- )
          _basisStatusRows[i] = rows[i];

        for( int j = numCols() - 1; j >= 0; j-- )
          _basisStatusCols[j] = cols[j];

        _hasBasis = true;
      }
  }



  /// clears starting basis
  template <class R>
  void SoPlexBase<R>::clearBasis()
  {
    _solver.reLoad();
    _status = _solver.status();
    _hasBasis = false;
    _rationalLUSolver.clear();
  }



  /// number of iterations since last call to solve
  template <class R>
	int SoPlexBase<R>::numIterations() const
  {
    return _statistics->iterations;
  }



  /// time spent in last call to solve
  template <class R>
  R SoPlexBase<R>::solveTime() const
  {
    return _statistics->solvingTime->time();
  }



  /// statistical information in form of a string
  template <class R>
  std::string SoPlexBase<R>::statisticString() const
  {
    std::stringstream s;
    s  << "Factorizations     : " << std::setw(10) << _statistics->luFactorizationsReal << std::endl
       << "  Time spent       : " << std::setw(10) << std::fixed << std::setprecision(2) << _statistics->luFactorizationTimeReal << std::endl
       << "Solves             : " << std::setw(10) << _statistics->luSolvesReal << std::endl
       << "  Time spent       : " << std::setw(10) << _statistics->luSolveTimeReal << std::endl
       << "Solution time      : " << std::setw(10) << std::fixed << std::setprecision(2) << solveTime() << std::endl
       << "Iterations         : " << std::setw(10) << numIterations() << std::endl;

    return s.str();
  }



  /// name of starter
  template <class R>
  const char* SoPlexBase<R>::getStarterName()
  {
    if( _starter )
      return _starter->getName();
    else
      return "none";
  }



  /// name of simplifier
  template <class R>
  const char* SoPlexBase<R>::getSimplifierName()
  {
    if( _simplifier )
      return _simplifier->getName();
    else
      return "none";
  }



  /// name of scaling method after simplifier
  template <class R>
  const char* SoPlexBase<R>::getScalerName()
  {
    if( _scaler )
      return _scaler->getName();
    else
      return "none";
  }



  /// name of currently loaded pricer
  template <class R>
  const char* SoPlexBase<R>::getPricerName()
  {
    return _solver.pricer()->getName();
  }



  /// name of currently loaded ratiotester
  template <class R>
  const char* SoPlexBase<R>::getRatiotesterName()
  {
    return _solver.ratiotester()->getName();
  }


  /// returns boolean parameter value
  template <class R>
	bool SoPlexBase<R>::boolParam(const BoolParam param) const
  {
    assert(param >= 0);
    assert(param < SoPlexBase<R>::BOOLPARAM_COUNT);
    return _currentSettings->_boolParamValues[param];
  }

  /// returns integer parameter value
  template <class R>
	int SoPlexBase<R>::intParam(const IntParam param) const
  {
    assert(param >= 0);
    assert(param < INTPARAM_COUNT);
    return _currentSettings->_intParamValues[param];
  }

  /// returns real parameter value
  template <class R>
  R SoPlexBase<R>::realParam(const RealParam param) const
  {
    assert(param >= 0);
    assert(param < REALPARAM_COUNT);
    return _currentSettings->_realParamValues[param];
  }

#ifdef SOPLEX_WITH_RATIONALPARAM
  /// returns rational parameter value
  Rational SoPlexBase<R>::rationalParam(const RationalParam param) const
  {
    assert(param >= 0);
    assert(param < RATIONALPARAM_COUNT);
    return _currentSettings->_rationalParamValues[param];
  }
#endif



  /// returns current parameter settings
  template <class R>
  const typename SoPlexBase<R>::Settings& SoPlexBase<R>::settings() const
  {
    return *_currentSettings;
  }



  /// sets boolean parameter value; returns true on success
  template <class R>
	bool SoPlexBase<R>::setBoolParam(const BoolParam param, const bool value, const bool init)
  {
    assert(param >= 0);
    assert(param < SoPlexBase<R>::BOOLPARAM_COUNT);
    assert(init || _isConsistent());

    if( !init && value == boolParam(param) )
      return true;

    switch( param )
      {
      case LIFTING:
        break;
      case EQTRANS:
        break;
      case TESTDUALINF:
        break;
      case RATFAC:
        break;
      case USEDECOMPDUALSIMPLEX:
        break;
      case COMPUTEDEGEN:
        break;
      case USECOMPDUAL:
        break;
      case EXPLICITVIOL:
        break;
      case ACCEPTCYCLING:
        break;
      case RATREC:
        break;
      case POWERSCALING:
        break;
      case RATFACJUMP:
        break;
      case ROWBOUNDFLIPS:
        _ratiotesterBoundFlipping.useBoundFlipsRow(value);
        break;
      case PERSISTENTSCALING:
        break;
      case FULLPERTURBATION:
        _solver.useFullPerturbation(value);
        break;
      case ENSURERAY:
        break;
      default:
        return false;
      }

    _currentSettings->_boolParamValues[param] = value;
    return true;
  }

  /// sets integer parameter value; returns true on success
  template <class R>
	bool SoPlexBase<R>::setIntParam(const IntParam param, const int value, const bool init)
  {
    assert(param >= 0);
    assert(param < INTPARAM_COUNT);
    assert(init || _isConsistent());

    if( !init && value == intParam(param) )
      return true;

    // check for a valid parameter value wrt bounds
    if( value < _currentSettings->intParam.lower[param] || value > _currentSettings->intParam.upper[param] )
      return false;

    switch( param )
      {
        // objective sense
      case SoPlexBase<R>::OBJSENSE:
        if( value != SoPlexBase<R>::OBJSENSE_MAXIMIZE && value != SoPlexBase<R>::OBJSENSE_MINIMIZE )
          return false;
        _realLP->changeSense(value == SoPlexBase<R>::OBJSENSE_MAXIMIZE ? SPxLPBase<R>::MAXIMIZE : SPxLPBase<R>::MINIMIZE);
        if( _rationalLP != 0 )
          _rationalLP->changeSense(value == SoPlexBase<R>::OBJSENSE_MAXIMIZE ? SPxLPRational::MAXIMIZE : SPxLPRational::MINIMIZE);
        _invalidateSolution();
        break;

        // type of computational form, i.e., column or row representation
      case SoPlexBase<R>::REPRESENTATION:
        if( value != SoPlexBase<R>::REPRESENTATION_COLUMN && value != SoPlexBase<R>::REPRESENTATION_ROW && value != SoPlexBase<R>::REPRESENTATION_AUTO )
          return false;
        break;

        // type of algorithm, i.e., primal or dual
      case SoPlexBase<R>::ALGORITHM:
        // decide upon entering/leaving at solve time depending on representation
        break;

        // type of LU update
      case SoPlexBase<R>::FACTOR_UPDATE_TYPE:
        if( value != SoPlexBase<R>::FACTOR_UPDATE_TYPE_ETA && value != SoPlexBase<R>::FACTOR_UPDATE_TYPE_FT )
          return false;
        _slufactor.setUtype(value == SoPlexBase<R>::FACTOR_UPDATE_TYPE_ETA ? SLUFactor::ETA : SLUFactor::FOREST_TOMLIN);
        break;

        // maximum number of updates before fresh factorization
      case SoPlexBase<R>::FACTOR_UPDATE_MAX:
        if( value == 0 )
          _solver.basis().setMaxUpdates(DEFAULT_REFACTOR_INTERVAL);
        else
          _solver.basis().setMaxUpdates(value);
        break;

        // iteration limit (-1 if unlimited)
      case SoPlexBase<R>::ITERLIMIT:
        break;

        // refinement limit (-1 if unlimited)
      case SoPlexBase<R>::REFLIMIT:
        break;

        // stalling refinement limit (-1 if unlimited)
      case SoPlexBase<R>::STALLREFLIMIT:
        break;

        // display frequency
      case SoPlexBase<R>::DISPLAYFREQ:
        _solver.setDisplayFreq(value);
        break;

        // verbosity level
      case SoPlexBase<R>::VERBOSITY:
        switch(value)
          {
          case 0:
            spxout.setVerbosity(SPxOut::ERROR);
            break;
          case 1:
            spxout.setVerbosity(SPxOut::WARNING);
            break;
          case 2:
            spxout.setVerbosity(SPxOut::DEBUG);
            break;
          case 3:
            spxout.setVerbosity(SPxOut::INFO1);
            break;
          case 4:
            spxout.setVerbosity(SPxOut::INFO2);
            break;
          case 5:
            spxout.setVerbosity(SPxOut::INFO3);
            break;
          }
        break;

        // type of simplifier
      case SoPlexBase<R>::SIMPLIFIER:
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
      case SoPlexBase<R>::SCALER:
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
          case SCALER_LEASTSQ:
            _scaler = &_scalerLeastsq;
            break;
          case SCALER_GEOEQUI:
            _scaler = &_scalerGeoequi;
            break;
          default:
            return false;
          }
        break;

        // type of starter used to create crash basis
      case SoPlexBase<R>::STARTER:
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
      case SoPlexBase<R>::PRICER:
        switch( value )
          {
          case PRICER_AUTO:
            _solver.setPricer(&_pricerAuto);
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
          default:
            return false;
          }
        break;

        // mode for synchronizing real and rational LP
      case SoPlexBase<R>::SYNCMODE:
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
      case SoPlexBase<R>::READMODE:
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
      case SoPlexBase<R>::SOLVEMODE:
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

        // mode for a posteriori feasibility checks; nothing to do but change the value if valid
      case SoPlexBase<R>::CHECKMODE:
        switch( value )
          {
          case CHECKMODE_REAL:
          case CHECKMODE_AUTO:
          case CHECKMODE_RATIONAL:
            break;
          default:
            return false;
          }
        break;

        // type of ratio test
      case SoPlexBase<R>::RATIOTESTER:
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

        // type of timer
      case SoPlexBase<R>::TIMER:
        switch( value )
          {
          case TIMER_OFF:
            _solver.setTiming( Timer::OFF);
            break;
          case TIMER_CPU:
            _solver.setTiming( Timer::USER_TIME );
            break;
          case TIMER_WALLCLOCK:
            _solver.setTiming( Timer::WALLCLOCK_TIME);
            break;
          default:
            return false;
          }
        break;

        // mode of hyper pricing
      case SoPlexBase<R>::HYPER_PRICING:
        switch( value )
          {
          case HYPER_PRICING_OFF:
          case HYPER_PRICING_AUTO:
          case HYPER_PRICING_ON:
            break;
          default:
            return false;
          }
        break;

        // minimum number of stalling refinements since last pivot to trigger rational factorization
      case SoPlexBase<R>::RATFAC_MINSTALLS:
        break;

        // maximum number of conjugate gradient iterations in least square scaling
      case SoPlexBase<R>::LEASTSQ_MAXROUNDS:
        if( _scaler )
          _scaler->setIntParam(value);
        break;

        // mode of solution polishing
      case SoPlexBase<R>::SOLUTION_POLISHING:
        switch( value )
          {
          case POLISHING_OFF:
            _solver.setSolutionPolishing(SPxSolverBase<R>::POLISH_OFF);
            break;
          case POLISHING_INTEGRALITY:
            _solver.setSolutionPolishing(SPxSolverBase<R>::POLISH_INTEGRALITY);
            break;
          case POLISHING_FRACTIONALITY:
            _solver.setSolutionPolishing(SPxSolverBase<R>::POLISH_FRACTIONALITY);
            break;
          default:
            return false;
          }
        break;

        // the decomposition based simplex parameter settings
      case DECOMP_ITERLIMIT:
        break;
      case DECOMP_MAXADDEDROWS:
        break;
      case DECOMP_DISPLAYFREQ:
        break;
      case DECOMP_VERBOSITY:
        break;

        // printing of condition n
      case PRINTBASISMETRIC:
        _solver.setMetricInformation  (value);
        break;

      default:
        return false;
      }

    _currentSettings->_intParamValues[param] = value;
    return true;
  }

  /// sets real parameter value; returns true on success
  template <class R>
	bool SoPlexBase<R>::setRealParam(const RealParam param, const R value, const bool init)
  {
    assert(param >= 0);
    assert(param < REALPARAM_COUNT);
    assert(init || _isConsistent());

    if( !init && value == realParam(param) )
      return true;

    if( value < _currentSettings->realParam.lower[param] || value > _currentSettings->realParam.upper[param] )
      return false;

    // required to set a different feastol or opttol
    Real tmp_value = value;

    switch( param )
      {
        // primal feasibility tolerance; passed to the floating point solver only when calling solve()
      case SoPlexBase<R>::FEASTOL:
#ifndef SOPLEX_WITH_GMP
        if( value < DEFAULT_EPS_PIVOT )
          {
            MSG_WARNING( spxout, spxout << "Cannot set feasibility tolerance to small value " << value << " without GMP - using " << DEFAULT_EPS_PIVOT << ".\n");
            tmp_value = DEFAULT_EPS_PIVOT;
            _rationalFeastol = DEFAULT_EPS_PIVOT;
            break;
          }
#endif
        _rationalFeastol = value;
        break;

        // dual feasibility tolerance; passed to the floating point solver only when calling solve()
      case SoPlexBase<R>::OPTTOL:
#ifndef SOPLEX_WITH_GMP
        if( value < DEFAULT_EPS_PIVOT )
          {
            MSG_WARNING( spxout, spxout << "Cannot set optimality tolerance to small value " << value << " without GMP - using " << DEFAULT_EPS_PIVOT << ".\n");
            tmp_value = DEFAULT_EPS_PIVOT;
            _rationalOpttol = DEFAULT_EPS_PIVOT;
            break;
          }
#endif
        _rationalOpttol = value;
        break;

        // general zero tolerance
      case SoPlexBase<R>::EPSILON_ZERO:
        Param::setEpsilon(value);
        break;

        // zero tolerance used in factorization
      case SoPlexBase<R>::EPSILON_FACTORIZATION:
        Param::setEpsilonFactorization(value);
        break;

        // zero tolerance used in update of the factorization
      case SoPlexBase<R>::EPSILON_UPDATE:
        Param::setEpsilonUpdate(value);
        break;

        // pivot zero tolerance used in factorization (declare numerical singularity for small LU pivots)
      case SoPlexBase<R>::EPSILON_PIVOT:
        Param::setEpsilonPivot(value);
        break;

        // infinity threshold
      case SoPlexBase<R>::INFTY:
        _rationalPosInfty = value;
        _rationalNegInfty = -value;
        if( intParam(SoPlexBase<R>::SYNCMODE) != SYNCMODE_ONLYREAL )
          _recomputeRangeTypesRational();
        break;

        // time limit in seconds (INFTY if unlimited)
      case SoPlexBase<R>::TIMELIMIT:
        break;

        // lower limit on objective value is set in solveReal()
      case SoPlexBase<R>::OBJLIMIT_LOWER:
        break;

        // upper limit on objective value is set in solveReal()
      case SoPlexBase<R>::OBJLIMIT_UPPER:
        break;

        // working tolerance for feasibility in floating-point solver
      case SoPlexBase<R>::FPFEASTOL:
        break;

        // working tolerance for optimality in floating-point solver
      case SoPlexBase<R>::FPOPTTOL:
        break;

        // maximum increase of scaling factors between refinements
      case SoPlexBase<R>::MAXSCALEINCR:
        _rationalMaxscaleincr = value;
        break;

        // lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)
      case SoPlexBase<R>::LIFTMINVAL:
        break;

        // upper threshold in lifting (nonzero matrix coefficients with larger absolute value will be reformulated)
      case SoPlexBase<R>::LIFTMAXVAL:
        break;

        // threshold for sparse pricing
      case SoPlexBase<R>::SPARSITY_THRESHOLD:
        break;

        // threshold on number of rows vs. number of columns for switching from column to row representations in auto mode
      case SoPlexBase<R>::REPRESENTATION_SWITCH:
        break;

        // geometric frequency at which to apply rational reconstruction
      case SoPlexBase<R>::RATREC_FREQ:
        break;

        // minimal reduction (sum of removed rows/cols) to continue simplification
      case SoPlexBase<R>::MINRED:
        break;

      case SoPlexBase<R>::REFAC_BASIS_NNZ:
        break;

      case SoPlexBase<R>::REFAC_UPDATE_FILL:
        break;

      case SoPlexBase<R>::REFAC_MEM_FACTOR:
        break;

        // accuracy of conjugate gradient method in least squares scaling (higher value leads to more iterations)
      case SoPlexBase<R>::LEASTSQ_ACRCY:
        if( _scaler )
          _scaler->setRealParam(value);
        break;

        // objective offset
      case SoPlexBase<R>::OBJ_OFFSET:
        if( _realLP )
          _realLP->changeObjOffset(value);
        if( _rationalLP )
          _rationalLP->changeObjOffset(value);
        break;

      default:
        return false;
      }

    _currentSettings->_realParamValues[param] = tmp_value;
    return true;
  }


#ifdef SOPLEX_WITH_RATIONALPARAM
  /// sets rational parameter value; returns true on success
  template <class R>
	bool SoPlexBase<R>::setRationalParam(const RationalParam param, const Rational value, const bool init)
  {
    assert(param >= 0);
    assert(param < RATIONALPARAM_COUNT);
    assert(init || _isConsistent());

    if( !init && value == rationalParam(param) )
      return true;

    if( value < _currentSettings->rationalParam.lower[param] || value > _currentSettings->rationalParam.upper[param] )
      return false;

    switch( param )
      {
      default:
        // currently, there are no rational-valued parameters
        return false;
      }

    _currentSettings->_rationalParamValues[param] = value;
    return true;
  }
#endif



  /// sets parameter settings; returns true on success
  template <class R>
	bool SoPlexBase<R>::setSettings(const Settings& newSettings, const bool init)
  {
    assert(init || _isConsistent());

    bool success = true;

    *_currentSettings = newSettings;

    for( int i = 0; i < SoPlexBase<R>::BOOLPARAM_COUNT; i++ )
      success &= setBoolParam((BoolParam)i, _currentSettings->_boolParamValues[i], init);

    for( int i = 0; i < SoPlexBase<R>::INTPARAM_COUNT; i++ )
      success &= setIntParam((IntParam)i, _currentSettings->_intParamValues[i], init);

    for( int i = 0; i < SoPlexBase<R>::REALPARAM_COUNT; i++ )
      success &= setRealParam((RealParam)i, _currentSettings->_realParamValues[i], init);

#ifdef SOPLEX_WITH_RATIONALPARAM
    for( int i = 0; i < SoPlexBase<R>::RATIONALPARAM_COUNT; i++ )
      success &= setRationalParam((RationalParam)i, _currentSettings->_rationalParamValues[i], init);
#endif

    assert(_isConsistent());

    return success;
  }

  /// resets default parameter settings
  template <class R>
  void SoPlexBase<R>::resetSettings(const bool quiet, const bool init)
  {
    for( int i = 0; i < SoPlexBase<R>::BOOLPARAM_COUNT; i++ )
      setBoolParam((BoolParam)i, _currentSettings->boolParam.defaultValue[i], init);

    for( int i = 0; i < SoPlexBase<R>::INTPARAM_COUNT; i++ )
      setIntParam((IntParam)i, _currentSettings->intParam.defaultValue[i], init);

    for( int i = 0; i < SoPlexBase<R>::REALPARAM_COUNT; i++ )
      setRealParam((RealParam)i, _currentSettings->realParam.defaultValue[i], init);

#ifdef SOPLEX_WITH_RATIONALPARAM
    for( int i = 0; i < SoPlexBase<R>::RATIONALPARAM_COUNT; i++ )
      success &= setRationalParam((RationalParam)i, _currentSettings->rationalParam.defaultValue[i], init);
#endif
  }

  /// print non-default parameter values
  template <class R>
  void SoPlexBase<R>::printUserSettings()
  {
    bool printedValue = false;

    SPxOut::setFixed(spxout.getCurrentStream());

    for( int i = 0; i < SoPlexBase<R>::BOOLPARAM_COUNT; i++ )
      {
        if( _currentSettings->_boolParamValues[i] == _currentSettings->boolParam.defaultValue[i] )
          continue;

        spxout << "bool:" << _currentSettings->boolParam.name[i] << " = " << (_currentSettings->_boolParamValues[i] ? "true\n" : "false\n");
        printedValue = true;
      }

    for( int i = 0; i < SoPlexBase<R>::INTPARAM_COUNT; i++ )
      {
        if( _currentSettings->_intParamValues[i] == _currentSettings->intParam.defaultValue[i] )
          continue;

        spxout << "int:" << _currentSettings->intParam.name[i] << " = " << _currentSettings->_intParamValues[i] << "\n";
        printedValue = true;
      }

    SPxOut::setScientific(spxout.getCurrentStream());

    for( int i = 0; i < SoPlexBase<R>::REALPARAM_COUNT; i++ )
      {
        if( _currentSettings->_realParamValues[i] == _currentSettings->realParam.defaultValue[i] )
          continue;

        spxout << "real:" << _currentSettings->realParam.name[i] << " = " << _currentSettings->_realParamValues[i] << "\n";
        printedValue = true;
      }

#ifdef SOPLEX_WITH_RATIONALPARAM
    for( int i = 0; i < SoPlexBase<R>::RATIONALPARAM_COUNT; i++ )
      {
        if( _currentSettings->_rationalParamValues[i] == _currentSettings->rationalParam.defaultValue[i] )
          continue;

        spxout << "rational:" << _currentSettings->rationalParam.name[i] << " = " << _currentSettings->_rationalParamValues[i] << "\n";
        printedValue = true;
      }
#endif

    if( _solver.random.getSeed() != DEFAULT_RANDOM_SEED )
      {
        spxout << "uint:random_seed = " << _solver.random.getSeed() << "\n";
        printedValue = true;
      }

    if( printedValue )
      spxout << std::endl;
  }

  /// prints status
  template <class R>
  void SoPlexBase<R>::printStatus(std::ostream& os, typename SPxSolverBase<R>::Status stat)
  {
    os << "SoPlex status       : ";

    switch( stat )
      {
      case SPxSolverBase<R>::ERROR:
        os << "error [unspecified]";
        break;
      case SPxSolverBase<R>::NO_RATIOTESTER:
        os << "error [no ratiotester loaded]";
        break;
      case SPxSolverBase<R>::NO_PRICER:
        os << "error [no pricer loaded]";
        break;
      case SPxSolverBase<R>::NO_SOLVER:
        os << "error [no linear solver loaded]";
        break;
      case SPxSolverBase<R>::NOT_INIT:
        os << "error [not initialized]";
        break;
      case SPxSolverBase<R>::ABORT_CYCLING:
        os << "solving aborted [cycling]";
        break;
      case SPxSolverBase<R>::ABORT_TIME:
        os << "solving aborted [time limit reached]";
        break;
      case SPxSolverBase<R>::ABORT_ITER:
        os << "solving aborted [iteration limit reached]";
        break;
      case SPxSolverBase<R>::ABORT_VALUE:
        os << "solving aborted [objective limit reached]";
        break;
      case SPxSolverBase<R>::NO_PROBLEM:
        os << "no problem loaded";
        break;
      case SPxSolverBase<R>::REGULAR:
        os << "basis is regular";
        break;
      case SPxSolverBase<R>::SINGULAR:
        os << "basis is singular";
        break;
      case SPxSolverBase<R>::OPTIMAL:
        os << "problem is solved [optimal]";
        break;
      case SPxSolverBase<R>::UNBOUNDED:
        os << "problem is solved [unbounded]";
        break;
      case SPxSolverBase<R>::INFEASIBLE:
        os << "problem is solved [infeasible]";
        break;
      case SPxSolverBase<R>::INForUNBD:
        os << "problem is solved [infeasible or unbounded]";
        break;
      case SPxSolverBase<R>::OPTIMAL_UNSCALED_VIOLATIONS:
         os << "problem is solved [optimal with unscaled violations]";
         break;
      default:
      case SPxSolverBase<R>::UNKNOWN:
        os << "unknown";
        break;
      }

    os << "\n";
  }


  /// prints version and compilation options
  template <class R>
  void SoPlexBase<R>::printVersion() const
  {
    // do not use preprocessor directives within the MSG_INFO1 macro
#if (SOPLEX_SUBVERSION > 0)
    MSG_INFO1( spxout, spxout << "SoPlex version " << SOPLEX_VERSION/100
               << "." << (SOPLEX_VERSION % 100)/10
               << "." << SOPLEX_VERSION % 10
               << "." << SOPLEX_SUBVERSION );
#else
    MSG_INFO1( spxout, spxout << "SoPlex version " << SOPLEX_VERSION/100
               << "." << (SOPLEX_VERSION % 100)/10
               << "." << SOPLEX_VERSION % 10 );
#endif

#ifndef NDEBUG
    MSG_INFO1( spxout, spxout << " [mode: debug]" );
#else
    MSG_INFO1( spxout, spxout << " [mode: optimized]" );
#endif

    MSG_INFO1( spxout, spxout << " [precision: " << (int)sizeof(Real) << " byte]" );

#ifdef SOPLEX_WITH_GMP
#ifdef mpir_version
    MSG_INFO1( spxout, spxout << " [rational: MPIR " << mpir_version << "]" );
#else
    MSG_INFO1( spxout, spxout << " [rational: GMP " << gmp_version << "]" );
#endif
#else
    MSG_INFO1( spxout, spxout << " [rational: long double]" );
#endif

    MSG_INFO1( spxout, spxout << " [githash: " << getGitHash() << "]\n" );
  }


  /// checks if real LP and rational LP are in sync; dimensions will always be compared,
  /// vector and matrix values only if the respective parameter is set to true.
  /// If quiet is set to true the function will only display which vectors are different.
  template <class R>
	bool SoPlexBase<R>::areLPsInSync(const bool checkVecVals, const bool checkMatVals, const bool quiet) const
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
        MSG_INFO1( spxout, spxout << "The number of Rows in the Real LP does not match the one in the Rational LP."
                   << " Real LP: " << _realLP->nRows() << "  Rational LP: " << _rationalLP->nRows() << std::endl);
        result = false;
        nRowsMatch = false;
      }

    // compare number of Columns
    if( _realLP->nCols() != _rationalLP->nCols() )
      {
        MSG_INFO1( spxout, spxout << "The number of Columns in the Real LP does not match the one in the Rational LP."
                   << " Real LP: " << _realLP->nCols() << "  Rational LP: " << _rationalLP->nCols() << std::endl);
        result = false;
        nColsMatch = false;
      }

    // compare number of nonZeros
    if( _realLP->nNzos() != _rationalLP->nNzos() )
      {
        MSG_INFO1( spxout, spxout << "The number of nonZeros in the Real LP does not match the one in the Rational LP."
                   << " Real LP: " << _realLP->nNzos() << "  Rational LP: " << _rationalLP->nNzos() << std::endl);
        result = false;
      }

    // compare the dimensions of the right hand side vectors
    if( _realLP->rhs().dim() != _rationalLP->rhs().dim() )
      {
        MSG_INFO1( spxout, spxout << "The dimension of the right hand side vector of the Real LP does not match the one of the Rational LP."
                   << " Real LP: " << _realLP->rhs().dim() << "  Rational LP: " << _rationalLP->rhs().dim() << std::endl);
        result = false;
        rhsDimMatch = false;

      }

    // compare the dimensions of the left hand side vectors
    if( _realLP->lhs().dim() != _rationalLP->lhs().dim() )
      {
        MSG_INFO1( spxout, spxout << "The dimension of the left hand side vector of the Real LP does not match the one of the Rational LP."
                   << " Real LP: " << _realLP->lhs().dim() << "  Rational LP: " << _rationalLP->lhs().dim() << std::endl);
        result = false;
        lhsDimMatch = false;
      }

    // compare the dimensions of the objective function vectors
    if( _realLP->maxObj().dim() != _rationalLP->maxObj().dim() )
      {
        MSG_INFO1( spxout, spxout << "The dimension of the objective function vector of the Real LP does not match the one of the Rational LP."
                   << " Real LP: " << _realLP->maxObj().dim() << "  Rational LP: " << _rationalLP->maxObj().dim() << std::endl);
        result = false;
        maxObjDimMatch = false;
      }

    // compare the sense
    if( (int)_realLP->spxSense() != (int)_rationalLP->spxSense() )
      {
        MSG_INFO1( spxout, spxout << "The objective function sense of the Real LP does not match the one of the Rational LP."
                   << " Real LP: " << (_realLP->spxSense() == SPxLPBase<R>::MINIMIZE ? "MIN" : "MAX")
                   << "  Rational LP: " << (_rationalLP->spxSense() == SPxLPRational::MINIMIZE ? "MIN" : "MAX") << std::endl);
        result = false;
      }

    // compare the dimensions of upper bound vectors
    if( _realLP->upper().dim() != _rationalLP->upper().dim() )
      {
        MSG_INFO1( spxout, spxout << "The dimension of the upper bound vector of the Real LP does not match the one of the Rational LP."
                   << " Real LP: " << _realLP->upper().dim() << "  Rational LP: " << _rationalLP->upper().dim() << std::endl);
        result = false;
        upperDimMatch = false;
      }

    // compare the dimensions of the objective function vectors
    if( _realLP->lower().dim() != _rationalLP->lower().dim() )
      {
        MSG_INFO1( spxout, spxout << "The dimension of the lower bound vector of the Real LP does not match the one of the Rational LP."
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
                if( (GE(_realLP->rhs()[i], realParam(SoPlexBase<R>::INFTY)) != (_rationalLP->rhs()[i] >= _rationalPosInfty))
                    || (LT(_realLP->rhs()[i], realParam(SoPlexBase<R>::INFTY)) && _rationalLP->rhs()[i] < _rationalPosInfty
                        && !_rationalLP->rhs()[i].isAdjacentTo((double)_realLP->rhs()[i])) )
                  {
                    if( !quiet )
                      {
                        MSG_INFO1( spxout, spxout << "Entries number " << i << " of the right hand side vectors don't match."
                                   << " Real LP: " << _realLP->rhs()[i] << "  Rational LP: " << _rationalLP->rhs()[i] << std::endl);
                      }
                    rhsValMatch = false;
                    result = false;
                  }
              }

            if( !rhsValMatch && quiet )
              {
                MSG_INFO1( spxout, spxout << "The values of the right hand side vectors don't match." << std::endl );
              }
          }

        // compares the values of the left hand side vectors
        if( lhsDimMatch )
          {
            for( int i = 0; i < _realLP->lhs().dim(); i++ )
              {
                if( (LE(_realLP->lhs()[i], -realParam(SoPlexBase<R>::INFTY)) != (_rationalLP->lhs()[i] <= _rationalNegInfty))
                    || (GT(_realLP->lhs()[i], -realParam(SoPlexBase<R>::INFTY)) && _rationalLP->lhs()[i] > _rationalNegInfty
                        && !_rationalLP->lhs()[i].isAdjacentTo((double)_realLP->lhs()[i])) )
                  {
                    if( !quiet )
                      {
                        MSG_INFO1( spxout, spxout << "Entries number " << i << " of the left hand side vectors don't match."
                                   << " Real LP: " << _realLP->lhs()[i] << "  Rational LP: " << _rationalLP->lhs()[i] << std::endl);
                      }
                    lhsValMatch = false;
                    result = false;
                  }
              }

            if( !lhsValMatch && quiet )
              {
                MSG_INFO1( spxout, spxout << "The values of the left hand side vectors don't match." << std::endl );
              }
          }

        // compares the values of the objective function vectors
        if( maxObjDimMatch )
          {
            for( int i = 0; i < _realLP->maxObj().dim(); i++ )
              {
                if( !_rationalLP->maxObj()[i].isAdjacentTo((double)_realLP->maxObj()[i]) )
                  {
                    if( !quiet )
                      {
                        MSG_INFO1( spxout, spxout << "Entries number " << i << " of the objective function vectors don't match."
                                   << " Real LP: " << _realLP->maxObj()[i] << "  Rational LP: " << _rationalLP->maxObj()[i] << std::endl);
                      }
                    maxObjValMatch = false;
                    result = false;
                  }
              }

            if( !maxObjValMatch && quiet )
              {
                MSG_INFO1( spxout, spxout << "The values of the objective function vectors don't match." << std::endl );
              }
          }

        // compares the values of the upper bound vectors
        if( upperDimMatch )
          {
            for( int i = 0; i < _realLP->upper().dim(); i++ )
              {
                if( (GE(_realLP->upper()[i], realParam(SoPlexBase<R>::INFTY)) != (_rationalLP->upper()[i] >= _rationalPosInfty))
                    || (LT(_realLP->upper()[i], realParam(SoPlexBase<R>::INFTY)) && _rationalLP->upper()[i] < _rationalPosInfty
                        && !_rationalLP->upper()[i].isAdjacentTo((double)_realLP->upper()[i])) )
                  {
                    if( !quiet )
                      {
                        MSG_INFO1( spxout, spxout << "Entries number " << i << " of the upper bound vectors don't match."
                                   << " Real LP: " << _realLP->upper()[i] << "  Rational LP: " << _rationalLP->upper()[i] << std::endl);
                      }
                    upperValMatch = false;
                    result = false;
                  }
              }

            if( !upperValMatch && quiet )
              {
                MSG_INFO1( spxout, spxout << "The values of the upper bound vectors don't match." << std::endl );
              }
          }

        // compares the values of the lower bound vectors
        if( lowerDimMatch )
          {
            for( int i = 0; i < _realLP->lower().dim(); i++ )
              {
                if( (LE(_realLP->lower()[i], -realParam(SoPlexBase<R>::INFTY)) != (_rationalLP->lower()[i] <= _rationalNegInfty))
                    || (GT(_realLP->lower()[i], -realParam(SoPlexBase<R>::INFTY)) && _rationalLP->lower()[i] > _rationalNegInfty
                        && !_rationalLP->lower()[i].isAdjacentTo((double)_realLP->lower()[i])) )
                  {
                    if( !quiet )
                      {
                        MSG_INFO1( spxout, spxout << "Entries number " << i << " of the lower bound vectors don't match."
                                   << " Real LP: " << _realLP->lower()[i] << "  Rational LP: " << _rationalLP->lower()[i] << std::endl);
                      }
                    lowerValMatch = false;
                    result = false;
                  }
              }

            if( !lowerValMatch && quiet )
              {
                MSG_INFO1( spxout, spxout << "The values of the lower bound vectors don't match." << std::endl );
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
                if( !_rationalLP->colVector(i)[j].isAdjacentTo((double)_realLP->colVector(i)[j]) )
                  {
                    if( !quiet )
                      {
                        MSG_INFO1( spxout, spxout << "Entries number " << j << " of column number " << i << " don't match."
                                   << " Real LP: " << _realLP->colVector(i)[j] << "  Rational LP: " << _rationalLP->colVector(i)[j] << std::endl);
                      }
                    matrixValMatch = false;
                    result = false;
                  }
              }
          }

        if( !matrixValMatch && quiet )
          {
            MSG_INFO1( spxout, spxout << "The values of the matrices don't match." << std::endl );
          }
      }

    return result;
  }



  /// set the random seed of the solver instance
  template <class R>
  void SoPlexBase<R>::setRandomSeed(unsigned int seed)
  {
    _solver.random.setSeed(seed);
  }



  /// returns the current random seed of the solver instance or the one stored in the settings
  template <class R>
  unsigned int  SoPlexBase<R>::randomSeed() const
  {
    return _solver.random.getSeed();
  }



  /// extends sparse vector to hold newmax entries if and only if it holds no more free entries
  template <class R>
  void SoPlexBase<R>::_ensureDSVectorRationalMemory(DSVectorRational& vec, const int newmax) const
  {
    assert(newmax > vec.size());
    if( vec.size() >= vec.max() )
      vec.setMax(newmax);
  }



  /// creates a permutation for removing rows/columns from an array of indices
  template <class R>
  void SoPlexBase<R>::_idxToPerm(int* idx, int idxSize, int* perm, int permSize) const
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
  template <class R>
  void SoPlexBase<R>::_rangeToPerm(int start, int end, int* perm, int permSize) const
  {
    assert(perm != 0);
    assert(permSize >= 0);

    for( int i = 0; i < permSize; i++ )
      perm[i] = (i < start || i > end) ? i : -1;
  }



  /// checks consistency
  template <class R>
	bool SoPlexBase<R>::_isConsistent() const
  {
    assert(_statistics != 0);
    assert(_currentSettings != 0);

    assert(_realLP != 0);
    assert(_rationalLP != 0 || intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_ONLYREAL);

    assert(_realLP != &_solver || _isRealLPLoaded);
    assert(_realLP == &_solver || !_isRealLPLoaded);

    assert(!_hasBasis || _isRealLPLoaded || _basisStatusRows.size() == numRows());
    assert(!_hasBasis || _isRealLPLoaded || _basisStatusCols.size() == numCols());
    assert(_rationalLUSolver.status() == SLinSolverRational::UNLOADED || _hasBasis);
    assert(_rationalLUSolver.status() == SLinSolverRational::UNLOADED || _rationalLUSolver.dim() == _rationalLUSolverBind.size());
    assert(_rationalLUSolver.status() == SLinSolverRational::UNLOADED || _rationalLUSolver.dim() == numRowsRational());

    assert(_rationalLP == 0 || _colTypes.size() == numColsRational());
    assert(_rationalLP == 0 || _rowTypes.size() == numRowsRational());

    return true;
  }



  /// should solving process be stopped?
  template <class R>
	bool SoPlexBase<R>::_isSolveStopped(bool& stoppedTime, bool& stoppedIter) const
  {
    assert(_statistics != 0);

    stoppedTime = (realParam(TIMELIMIT) < realParam(INFTY) && _statistics->solvingTime->time() >= realParam(TIMELIMIT));
    stoppedIter = (intParam(ITERLIMIT) >= 0 && _statistics->iterations >= intParam(ITERLIMIT))
      || (intParam(REFLIMIT) >= 0 && _statistics->refinements >= intParam(REFLIMIT))
      || (intParam(STALLREFLIMIT) >= 0 && _statistics->stallRefinements >= intParam(STALLREFLIMIT));

    return stoppedTime || stoppedIter;
  }



  /// determines RangeType from real bounds
  template <class R>
  typename SoPlexBase<R>::RangeType SoPlexBase<R>::_rangeTypeReal(const R& lower, const R& upper) const
  {
    assert(lower <= upper);

    if( lower <= -infinity )
      {
        if( upper >= infinity )
          return RANGETYPE_FREE;
        else
          return RANGETYPE_UPPER;
      }
    else
      {
        if( upper >= infinity )
          return RANGETYPE_LOWER;
        else if( lower == upper )
          return RANGETYPE_FIXED;
        else
          return RANGETYPE_BOXED;
      }
  }


