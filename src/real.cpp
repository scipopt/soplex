/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2001 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: real.cpp,v 1.3 2002/01/30 16:22:18 bzfkocht Exp $"

/**@file  real.cpp
 * @brief Floating point type definition.
 */
#include "real.h"

namespace soplex
{

const Real infinity = DEFAULT_INFINITY;

Real Param::s_epsilon = DEFAULT_EPS_ZERO;

void Param::setEpsilon(Real eps)
{
   s_epsilon = eps;
}

#if 0
// This results (correctly) in a exception on alpha processors
void Param::computeEpsilon()
{
   volatile Real one = 1.0;
   volatile Real x;
   volatile Real store;

   for(x = one; store != one; x /= 10.0)
      store = one + x;

   s_epsilon = x / 100.0;
}
#endif

} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------

