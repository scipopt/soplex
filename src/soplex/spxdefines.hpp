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

template <class R>
void Param<R>::setEpsilonFactorization(R eps)
{
   s_epsilon_factorization = eps;
}


/// returns \c true iff |a-b| <= eps
  template <class R>
  inline bool EQ(Real a, R b, R eps = Param<R>::epsilon())
{
   return spxAbs(a - b) <= eps;
}

/// returns \c true iff |a-b| > eps
template <class R>
inline bool NE(Real a, R b, R eps = Param<R>::epsilon())
{
   return spxAbs(a - b) > eps;
}

/// returns \c true iff a < b + eps
template <class R>
inline bool LT(Real a, R b, R eps = Param<R>::epsilon())
{
   return (a - b) < -eps;
}

/// returns \c true iff a <= b + eps
template <class R>
inline bool LE(Real a, R b, R eps = Param<R>::epsilon())
{
   return (a - b) < eps;
}

/// returns \c true iff a > b + eps
template <class R>
inline bool GT(Real a, R b, R eps = Param<R>::epsilon())
{
   return (a - b) > eps;
}

/// returns \c true iff a >= b + eps
template <class R>
inline bool GE(Real a, R b, R eps = Param<R>::epsilon())
{
   return (a - b) > -eps;
}

/// returns \c true iff |a| <= eps
template <class R>
inline bool isZero(Real a, R eps = Param<R>::epsilon())
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
inline bool EQrel(Real a, R b, R eps = Param<R>::epsilon())
{
   return spxAbs(relDiff(a, b)) <= eps;
}

/// returns \c true iff |relDiff(a,b)| > eps
template <class R>
inline bool NErel(Real a, R b, R eps = Param<R>::epsilon())
{
   return spxAbs(relDiff(a, b)) > eps;
}

/// returns \c true iff relDiff(a,b) <= -eps
template <class R>
inline bool LTrel(Real a, R b, R eps = Param<R>::epsilon())
{
   return relDiff(a, b) <= -eps;
}

/// returns \c true iff relDiff(a,b) <= eps
template <class R>
inline bool LErel(Real a, R b, R eps = Param<R>::epsilon())
{
   return relDiff(a, b) <= eps;
}

/// returns \c true iff relDiff(a,b) > eps
template <class R>
inline bool GTrel(Real a, R b, R eps = Param<R>::epsilon())
{
   return relDiff(a, b) > eps;
}

/// returns \c true iff relDiff(a,b) > -eps
template <class R>
inline bool GErel(Real a, R b, R eps = Param<R>::epsilon())
{
   return relDiff(a, b) > -eps;
}

