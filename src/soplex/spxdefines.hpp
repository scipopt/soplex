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

/**@file  spxdefines.hpp
 * @brief General templated functions for SoPlex
 */

// Defining the static members of the Param class
// THREADLOCAL is a #define to thread_local. (Is it really needed?)
template <class R>
THREADLOCAL R Param<R>::s_epsilon               = DEFAULT_EPS_ZERO;
template <class R>
THREADLOCAL R Param<R>::s_epsilon_factorization = DEFAULT_EPS_FACTOR;
template <class R>
THREADLOCAL R Param<R>::s_epsilon_update        = DEFAULT_EPS_UPDATE;
template <class R>
THREADLOCAL R Param<R>::s_epsilon_pivot         = DEFAULT_EPS_PIVOT;


template <class R>
void Param<R>::setEpsilonFactorization(R eps)
{
   s_epsilon_factorization = eps;
}


/// returns \c true iff |a-b| <= eps
  template <class R>
  inline bool EQ(R a, R b, R eps = Param<R>::epsilon())
{
   return spxAbs(a - b) <= eps;
}

/// returns \c true iff |a-b| > eps
template <class R>
inline bool NE(R a, R b, R eps = Param<R>::epsilon())
{
   return spxAbs(a - b) > eps;
}

/// returns \c true iff a < b + eps
template <class R>
inline bool LT(R a, R b, R eps = Param<R>::epsilon())
{
   return (a - b) < -eps;
}

/// returns \c true iff a <= b + eps
template <class R>
inline bool LE(R a, R b, R eps = Param<R>::epsilon())
{
   return (a - b) < eps;
}

/// returns \c true iff a > b + eps
template <class R>
inline bool GT(R a, R b, R eps = Param<R>::epsilon())
{
   return (a - b) > eps;
}

/// returns \c true iff a >= b + eps
template <class R>
inline bool GE(R a, R b, R eps = Param<R>::epsilon())
{
   return (a - b) > -eps;
}

/// returns \c true iff |a| <= eps
template <class R>
inline bool isZero(R a, R eps = Param<R>::epsilon())
{
   return spxAbs(a) <= eps;
}

/// returns \c true iff |a| > eps
template <class R>
inline bool isNotZero(R a, R eps = Param<R>::epsilon())
{
   return spxAbs(a) > eps;
}

/// returns \c true iff |relDiff(a,b)| <= eps
template <class R>
inline bool EQrel(R a, R b, R eps = Param<R>::epsilon())
{
   return spxAbs(relDiff(a, b)) <= eps;
}

/// returns \c true iff |relDiff(a,b)| > eps
template <class R>
inline bool NErel(R a, R b, R eps = Param<R>::epsilon())
{
   return spxAbs(relDiff(a, b)) > eps;
}

/// returns \c true iff relDiff(a,b) <= -eps
template <class R>
inline bool LTrel(R a, R b, R eps = Param<R>::epsilon())
{
   return relDiff(a, b) <= -eps;
}

/// returns \c true iff relDiff(a,b) <= eps
template <class R>
inline bool LErel(R a, R b, R eps = Param<R>::epsilon())
{
   return relDiff(a, b) <= eps;
}

/// returns \c true iff relDiff(a,b) > eps
template <class R>
inline bool GTrel(R a, R b, R eps = Param<R>::epsilon())
{
   return relDiff(a, b) > eps;
}

/// returns \c true iff relDiff(a,b) > -eps
template <class R>
inline bool GErel(R a, R b, R eps = Param<R>::epsilon())
{
   return relDiff(a, b) > -eps;
}

// Templated functions, originally from spxdefines.cpp
template <class R>
 R Param<R>::epsilon()
{
  return R(s_epsilon);
}

template <class R>
void Param<R>::setEpsilon(R eps)
{
   s_epsilon = eps;
}

template <class R>
 R Param<R>::epsilonFactorization()
{
   return s_epsilon_factorization;
}

template <class R>
 R Param<R>::epsilonUpdate()
{
   return s_epsilon_update;
}

template <class R>
void Param<R>::setEpsilonUpdate(R eps)
{
  s_epsilon_update = R(eps);
}

template <class R>
 R Param<R>::epsilonPivot()
{
   return s_epsilon_pivot;
}

template <class R>
void Param<R>::setEpsilonPivot(R eps)
{
   s_epsilon_pivot = eps;
}
