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

/**@file  spxlpbase_real.hpp
 * @brief General templated function definitions for spxlpbase real
 */


/// Returns unscaled lower bound of column \p i.
template<class R>
R SPxLPBase<R>::lowerUnscaled(int i) const
{
  assert(i >= 0 && i < nCols());
  if( _isScaled )
    return lp_scaler->lowerUnscaled(*this, i);
  else
    return LPColSetBase<R>::lower(i);
}

