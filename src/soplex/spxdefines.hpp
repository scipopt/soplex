/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  spxdefines.hpp
 * @brief General templated functions for SoPlex
 */

// Defining the static members of the Param class
// THREADLOCAL is a #define to thread_local. (Is it really needed?)

namespace soplex
{

//   template <class R>
// THREADLOCAL R Param::s_epsilon               = DEFAULT_EPS_ZERO;

//   template <class R>
// THREADLOCAL R Param::s_epsilon_factorization = DEFAULT_EPS_FACTOR;

//   template <class R>
// THREADLOCAL R Param::s_epsilon_update        = DEFAULT_EPS_UPDATE;

//   template <class R>
// THREADLOCAL R Param::s_epsilon_pivot         = DEFAULT_EPS_PIVOT;


/// returns \c true iff |a-b| <= eps
template <class R, class S>
inline bool EQ(R a, S b, R eps = Param::epsilon())
{
   return spxAbs(R(a - b)) <= eps;
}

/// returns \c true iff |a-b| > eps
template <class R, class S>
inline bool NE(R a, S b, R eps = Param::epsilon())
{
   return spxAbs(a - b) > eps;
}

/// returns \c true iff a < b + eps
template <class R, class S>
inline bool LT(R a, S b, R eps = Param::epsilon())
{
   return (a - b) < -eps;
}

/// returns \c true iff a <= b + eps
template <class R, class S>
inline bool LE(R a, S b, R eps = Param::epsilon())
{
   return (a - b) < eps;
}

/// returns \c true iff a > b + eps
template <class R, class S>
inline bool GT(R a, S b, R eps = Param::epsilon())
{
   return (a - b) > eps;
}

/// returns \c true iff a >= b + eps
template <class R, class S>
inline bool GE(R a, S b, R eps = Param::epsilon())
{
   return (a - b) > -eps;
}

/// returns \c true iff |a| <= eps
template <class R>
inline bool isZero(R a, R eps = Param::epsilon())
{
   return spxAbs(a) <= eps;
}

/// returns \c true iff |a| > eps
template <class R>
inline bool isNotZero(R a, R eps = Param::epsilon())
{
   return spxAbs(a) > eps;
}

/// returns \c true iff |relDiff(a,b)| <= eps
template <class R, class S>
inline bool EQrel(R a, S b, R eps = Param::epsilon())
{
   return spxAbs(relDiff(a, b)) <= eps;
}

/// returns \c true iff |relDiff(a,b)| > eps
template <class R, class S>
inline bool NErel(R a, S b, R eps = Param::epsilon())
{
   return spxAbs(relDiff(a, b)) > eps;
}

/// returns \c true iff relDiff(a,b) <= -eps
template <class R, class S>
inline bool LTrel(R a, S b, R eps = Param::epsilon())
{
   return relDiff(a, b) <= -eps;
}

/// returns \c true iff relDiff(a,b) <= eps
template <class R, class S>
inline bool LErel(R a, S b, R eps = Param::epsilon())
{
   return relDiff(a, b) <= eps;
}

/// returns \c true iff relDiff(a,b) > eps
template <class R, class S>
inline bool GTrel(R a, S b, R eps = Param::epsilon())
{
   return relDiff(a, b) > eps;
}

/// returns \c true iff relDiff(a,b) > -eps
template <class R, class S>
inline bool GErel(R a, S b, R eps = Param::epsilon())
{
   return relDiff(a, b) > -eps;
}

// Instantiation for Real
inline Real spxLdexp(Real x, int exp)
{
   return std::ldexp(x, exp);
}

inline Real spxFrexp(Real x, int* exp)
{
   return std::frexp(x, exp);
}

} // namespace soplex
