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
#pragma ident "@(#) $Id: real.h,v 1.6 2002/01/29 15:38:48 bzfkocht Exp $"

/**@file  real.h
 * @brief Floating point type definition.
 */
#ifndef _REAL_H_
#define _REAL_H_

#include <math.h>

namespace soplex
{
#ifdef WITH_LONG_DOUBLE

typedef long double Real;

#define DEFAULT_EPS_ZERO  1e-32  // additive zero. 1.0 + EPS_ZERO == 1.0
#define DEFAULT_INFINITY  1e100

#else

typedef double Real;

#define DEFAULT_EPS_ZERO  1e-16  // additive zero. 1.0 + EPS_ZERO == 1.0
#define DEFAULT_INFINITY  1e100

#endif // !WITH_LONG_DOUBLE

extern const Real eps_zero;
extern const Real infinity;

inline bool EQ(Real a, Real b, Real eps = eps_zero)
{
   return fabs(a - b) <= eps;
}

inline bool NE(Real a, Real b, Real eps = eps_zero)
{
   return fabs(a - b) > eps;
}

inline bool LT(Real a, Real b, Real eps = eps_zero)
{
   return (a - b) < -eps;
}

inline bool LE(Real a, Real b, Real eps = eps_zero)
{
   return (a - b) < eps;
}

inline bool GT(Real a, Real b, Real eps = eps_zero)
{
   return (a - b) > eps;
}

inline bool GE(Real a, Real b, Real eps = eps_zero)
{
   return (a - b) > -eps;
}

} // namespace soplex
#endif // _REAL_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------



