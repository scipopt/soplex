/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxdefines.h,v 1.20 2005/01/03 14:46:48 bzfkocht Exp $"

/**@file  spxdefines.h
 * @brief Debugging, floating point type and parameter definitions.
 *
 * In optimized code with \c NDEBUG defined, only #VERBOSE1, #VERBOSE2
 * and #VERBOSE3 are set.
 * If \c NDEBUG is not defined, the code within #TRACE is used.
 * If \c DEBUGGING is defined, the code within #DEBUG is also used.
 * If \c TRACE_METHOD is defined, the method tracing with #METHOD is
 * activated.
 *
 * If \c WITH_LONG_DOUBLE is defined, all #Real numbers are of type 
 * long double instead of just double.
 */
#ifndef _SPXDEFINES_H_
#define _SPXDEFINES_H_

#include <math.h>

#ifdef TRACE_METHOD
#include "tracemethod.h"
#endif


namespace soplex
{
#define SOPLEX_VERSION   122

/*-----------------------------------------------------------------------------
 * Debugging Macros etc.
 *-----------------------------------------------------------------------------
 */
#define VERBOSE1(x) { if(Param::verbose() >= 1) {x} }
#define VERBOSE2(x) { if(Param::verbose() >= 2) {x} }
#define VERBOSE3(x) { if(Param::verbose() >= 3) {x} }

#ifndef NDEBUG
#define TRACE(x) {x}
#else
#define TRACE(x) /**/
#endif //!NDEBUG

#if defined(DEBUGGING)
#define DEBUG(x) {x}
#else
#define DEBUG(x) /**/
#endif //!DEBUGGING

#if defined(TRACE_METHOD)
#define METHOD(x) soplex::TraceMethod _trace_method_(x, __FILE__, __LINE__)
#else
#define METHOD(x) /**/
#endif // !TRACE_METHOD

/*-----------------------------------------------------------------------------
 * Long double support, Parameters and Epsilons
 *-----------------------------------------------------------------------------
 */
#ifdef WITH_LONG_DOUBLE

typedef long double Real;

#ifndef REAL
#define REAL(x)  x##L
#endif
#ifndef DEFAULT_BND_VIOL
#define DEFAULT_BND_VIOL   1e-12
#endif
#ifndef DEFAULT_EPS_ZERO
#define DEFAULT_EPS_ZERO   1e-28  // ~ additive zero. 1.0 + EPS_ZERO == 1.0
#endif
#ifndef DEFAULT_EPS_FACTOR
#define DEFAULT_EPS_FACTOR 1e-30
#endif
#ifndef DEFAULT_EPS_UPDATE
#define DEFAULT_EPS_UPDATE 1e-26
#endif
#define DEFAULT_INFINITY   1e100

#else

typedef double Real;

#ifndef REAL
#define REAL(x)  x
#endif
#ifndef DEFAULT_BND_VIOL
#define DEFAULT_BND_VIOL   1e-6
#endif
#ifndef DEFAULT_EPS_ZERO
#define DEFAULT_EPS_ZERO   1e-16  // ~ additive zero. 1.0 + EPS_ZERO == 1.0
#endif
#ifndef DEFAULT_EPS_FACTOR
#define DEFAULT_EPS_FACTOR 1e-20
#endif
#ifndef DEFAULT_EPS_UPDATE
#define DEFAULT_EPS_UPDATE 1e-16
#endif
#define DEFAULT_INFINITY   1e100

#endif // !WITH_LONG_DOUBLE

extern const Real infinity;

class Param
{
private:
   static Real s_epsilon;
   static Real s_epsilon_factorization;
   static Real s_epsilon_update;
   static int  s_verbose;

public:
   ///
   inline static Real epsilon() 
   {
      return s_epsilon;
   }
   ///
   static void setEpsilon(Real eps);
   ///
   inline static Real epsilonFactorization()
   {
      return s_epsilon_factorization;
   }
   ///
   static void setEpsilonFactorization(Real eps);
   ///
   inline static Real epsilonUpdate()
   {
      return s_epsilon_update;
   }
   ///
   static void setEpsilonUpdate(Real eps);
   ///
   inline static int verbose()
   {
      return s_verbose;
   }
   static void setVerbose(int p_verbose);
};

// (max(|a|,|b|)
inline static Real maxAbs(Real a, Real b)
{
   return fabs(a) > fabs(b) ? fabs(a) : fabs(b);
}

// (a-b)/max(|a|,|b|,1.0)
inline static Real relDiff(Real a, Real b)
{
   Real scale = 1.0 / maxAbs(maxAbs(a, b), 1.0);

   return a * scale - b * scale;
}

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

inline bool EQrel(Real a, Real b, Real eps = Param::epsilon())
{
   return fabs(relDiff(a, b)) <= eps;
}

inline bool NErel(Real a, Real b, Real eps = Param::epsilon())
{
   return fabs(relDiff(a, b)) > eps;
}

inline bool LTrel(Real a, Real b, Real eps = Param::epsilon())
{
   return relDiff(a, b) < -eps;
}

inline bool LErel(Real a, Real b, Real eps = Param::epsilon())
{
   return relDiff(a, b) < eps;
}

inline bool GTrel(Real a, Real b, Real eps = Param::epsilon())
{
   return relDiff(a, b) > eps;
}

inline bool GErel(Real a, Real b, Real eps = Param::epsilon())
{
   return relDiff(a, b) > -eps;
}

} // namespace soplex
#endif // _SPXDEFINES_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------



