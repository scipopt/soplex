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
#pragma ident "@(#) $Id: real.h,v 1.7 2002/01/30 14:14:00 bzfkocht Exp $"

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

#define DEFAULT_BND_VIOL  1e-12
#define DEFAULT_EPS_ZERO  1e-32  // ~ additive zero. 1.0 + EPS_ZERO == 1.0
#define DEFAULT_INFINITY  1e100

#else

typedef double Real;

#define DEFAULT_BND_VIOL  1e-6
#define DEFAULT_EPS_ZERO  1e-18  // ~ additive zero. 1.0 + EPS_ZERO == 1.0
#define DEFAULT_INFINITY  1e100

#endif // !WITH_LONG_DOUBLE

extern const Real infinity;

class Param
{
private:
   static Real s_epsilon;

public:
   inline static Real epsilon()
   {
      return s_epsilon;
   }
   static void setEpsilon(Real eps);
   static void computeEpsilon();
};

inline bool EQ(Real a, Real b, Real eps = Param::epsilon())
{
   return fabs(a - b) <= eps;
}

inline bool NE(Real a, Real b, Real eps = Param::epsilon())
{
   return fabs(a - b) > eps;
}

inline bool LT(Real a, Real b, Real eps = Param::epsilon())
{
   return (a - b) < -eps;
}

inline bool LE(Real a, Real b, Real eps = Param::epsilon())
{
   return (a - b) < eps;
}

inline bool GT(Real a, Real b, Real eps = Param::epsilon())
{
   return (a - b) > eps;
}

inline bool GE(Real a, Real b, Real eps = Param::epsilon())
{
   return (a - b) > -eps;
}

inline bool isZero(Real a, Real eps = Param::epsilon())
{
   return fabs(a) <= eps;
}

inline bool isNotZero(Real a, Real eps = Param::epsilon())
{
   return fabs(a) > eps;
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



