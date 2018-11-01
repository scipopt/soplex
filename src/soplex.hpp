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

