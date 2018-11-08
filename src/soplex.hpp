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
Real SoPlexBase<R>::minAbsNonzeroReal() const
{
  assert(_realLP != 0);
  return _realLP->minAbsNzo();
}


/// returns biggest non-zero element in absolute value
template <class R>
Real SoPlexBase<R>::maxAbsNonzeroReal() const
{
  assert(_realLP != 0);
  return _realLP->maxAbsNzo();
}


/// returns (unscaled) coefficient
template <class R>
Real SoPlexBase<R>::coefReal(int row, int col) const
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
Real SoPlexBase<R>::rhsReal(int i) const
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
Real SoPlexBase<R>::lhsReal(int i) const
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
