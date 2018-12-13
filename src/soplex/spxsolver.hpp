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

