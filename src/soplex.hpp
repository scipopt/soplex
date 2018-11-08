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

