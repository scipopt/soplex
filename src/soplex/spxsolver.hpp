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

/**@file  spxsolver.hpp
 * @brief General templated functions for SoPlex
 */

template <class R>
bool SPxSolverBase<R>::read(std::istream& in, NameSet* rowNames,
                            NameSet* colNames, DIdxSet* intVars)
{
  if( initialized )
    {
      clear();
      unInit();

      if (thepricer)
        thepricer->clear();

      if (theratiotester)
        theratiotester->clear();
    }

  this->unLoad();
  if (!SPxLPBase<R>::read(in, rowNames, colNames, intVars))
    return false;

  this->theLP = this;

  return true;
}


template <class R>
void SPxSolverBase<R>::setFeastol(R d)
{

  if( d <= 0.0 )
    throw SPxInterfaceException("XSOLVE30 Cannot set feastol less than or equal to zero.");

  if( theRep == COLUMN )
    m_entertol = d;
  else
    m_leavetol = d;
}


template <class R>
void SPxSolverBase<R>::setDelta(R d)
{

  if( d <= 0.0 )
    throw SPxInterfaceException("XSOLVE32 Cannot set delta less than or equal to zero.");

  m_entertol = d;
  m_leavetol = d;
}

template <class R>
void SPxSolverBase<R>::loadLP(const SPxLPBase<R>& lp, bool initSlackBasis)
{
  clear();
  unInit();
  this->unLoad();
  resetClockStats();
  if (thepricer)
    thepricer->clear();
  if (theratiotester)
    theratiotester->clear();
  SPxLPBase<R>::operator=(lp);
  reDim();
  SPxBasisBase<R>::load(this, initSlackBasis);
}

template <class R>
void SPxSolverBase<R>::setTerminationValue(R p_value)
{
  objLimit = p_value;
}

template <class R>
void SPxSolverBase<R>::changeLower(int i, const R& newLower, bool scale)
{
  if( newLower != (scale ? SPxLPBase<R>::lowerUnscaled(i) : SPxLPBase<R>::lower(i)) )
    {
      R oldLower = SPxLPBase<R>::lower(i);
      // This has to be done before calling changeLowerStatus() because that is calling
      // basis.dualColStatus() which calls lower() and needs the changed value.
      SPxLPBase<R>::changeLower(i, newLower, scale);

      if (SPxBasisBase<R>::status() > SPxBasisBase<R>::NO_PROBLEM)
        {
          changeLowerStatus(i, SPxLPBase<R>::lower(i), oldLower);
          unInit();
        }
    }
}
template <class R>
SPxSolverBase<R>::SPxSolverBase(
                                Type            p_type,
                                Representation  p_rep,
                                Timer::TYPE     ttype)
  : theType (p_type)
  , thePricing(FULL)
  , theRep(p_rep)
  , polishObj(POLISH_OFF)
  , theTime(0)
  , timerType(ttype)
  , theCumulativeTime(0.0)
  , maxIters (-1)
  , maxTime (infinity)
  , nClckSkipsLeft(0)
  , nCallsToTimelim(0)
  , objLimit(infinity)
  , m_status(UNKNOWN)
  , m_nonbasicValue(0.0)
  , m_nonbasicValueUpToDate(false)
  , m_pricingViol(0.0)
  , m_pricingViolUpToDate(false)
  , m_pricingViolCo(0.0)
  , m_pricingViolCoUpToDate(false)
  , theShift (0)
  , m_maxCycle(100)
  , m_numCycle(0)
  , initialized (false)
  , solveVector2 (0)
  , solveVector3 (0)
  , coSolveVector2(0)
  , coSolveVector3(0)
  , freePricer (false)
  , freeRatioTester (false)
  , freeStarter (false)
  , displayLine (0)
  , displayFreq (200)
  , sparsePricingFactor(SPARSITYFACTOR)
  , getStartingDecompBasis(false)
  , computeDegeneracy(false)
  , degenCompIterOffset(0)
  , fullPerturbation(false)
  , printBasisMetric(0)
  , unitVecs (0)
  , primVec (0, Param<R>::epsilon())
  , dualVec (0, Param<R>::epsilon())
  , addVec (0, Param<R>::epsilon())
  , thepricer (0)
  , theratiotester (0)
  , thestarter (0)
  , boundrange(0.0)
  , siderange(0.0)
  , objrange(0.0)
  , infeasibilities(0)
  , infeasibilitiesCo(0)
  , isInfeasible(0)
  , isInfeasibleCo(0)
  , sparsePricingLeave(false)
  , sparsePricingEnter(false)
  , sparsePricingEnterCo(false)
  , hyperPricingLeave(true)
  , hyperPricingEnter(true)
  , remainingRoundsLeave(0)
  , remainingRoundsEnter(0)
  , remainingRoundsEnterCo(0)
  , weights(0)
  , coWeights(0)
  , weightsAreSetup(false)
  , integerVariables(0)
{
  theTime = TimerFactory::createTimer(timerType);

  setDelta(DEFAULT_BND_VIOL);

  this->theLP = this;
  initRep(p_rep);

  // info: SPxBasisBase is not consistent in this moment.
  //assert(SPxSolverBase<R>::isConsistent());
}

template <class R>
SPxSolverBase<R>::~SPxSolverBase()
{
  assert(!freePricer || thepricer != 0);
  assert(!freeRatioTester || theratiotester != 0);
  assert(!freeStarter || thestarter != 0);

  if(freePricer)
    {
      delete thepricer;
      thepricer = 0;
    }

  if(freeRatioTester)
    {
      delete theratiotester;
      theratiotester = 0;
    }

  if(freeStarter)
    {
      delete thestarter;
      thestarter = 0;
    }

  // free timer
  assert(theTime);
  theTime->~Timer();
  spx_free(theTime);
}

