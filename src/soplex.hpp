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

