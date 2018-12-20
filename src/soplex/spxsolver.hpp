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
