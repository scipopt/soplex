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
bool SoPlexBase<R>::getDualViolation(Real& maxviol, Real& sumviol)
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
bool SoPlexBase<R>::getRedCostViolation(Real& maxviol, Real& sumviol)
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
