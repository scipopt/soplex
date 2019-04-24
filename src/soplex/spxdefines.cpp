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

/**@file  spxdefines.cpp
 * @brief Debugging, floating point type and parameter definitions.
 */
#include "assert.h"
#include "soplex/spxdefines.h"
#include "soplex/spxout.h"
#include "soplex/rational.h"

#include "boost/multiprecision/number.hpp"
#include "boost/multiprecision/debug_adaptor.hpp"

namespace soplex
{
  // Overloaded EQ function
  bool EQ(int a, int b)
  {
    return (a == b);
  }

  THREADLOCAL const Real infinity                 = DEFAULT_INFINITY;

  THREADLOCAL Real Param::s_epsilon               = DEFAULT_EPS_ZERO;

  THREADLOCAL Real Param::s_epsilon_factorization = DEFAULT_EPS_FACTOR;

  THREADLOCAL Real Param::s_epsilon_update        = DEFAULT_EPS_UPDATE;

  THREADLOCAL Real Param::s_epsilon_pivot         = DEFAULT_EPS_PIVOT;

  using mpfr_float = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<50>, boost::multiprecision::et_off>;
  using mpfr_debug = boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::mpfr_float_backend<50> >, boost::multiprecision::et_off>;

  // // @todo: have this in the hpp file?
  // template <>
  // THREADLOCAL  mpfr_float Param<mpfr_float >::s_epsilon               = DEFAULT_EPS_ZERO;

  // template <>
  // THREADLOCAL mpfr_float Param<mpfr_float >::s_epsilon_factorization = DEFAULT_EPS_FACTOR;

  // template <>
  // THREADLOCAL mpfr_float Param<mpfr_float >::s_epsilon_update        = DEFAULT_EPS_UPDATE;

  // template <>
  // THREADLOCAL mpfr_float Param<mpfr_float >::s_epsilon_pivot         = DEFAULT_EPS_PIVOT;


  // template <>
  // THREADLOCAL  mpfr_debug Param<mpfr_debug >::s_epsilon               = DEFAULT_EPS_ZERO;

  // template <>
  // THREADLOCAL mpfr_debug Param<mpfr_debug >::s_epsilon_factorization = DEFAULT_EPS_FACTOR;

  // template <>
  // THREADLOCAL mpfr_debug Param<mpfr_debug >::s_epsilon_update        = DEFAULT_EPS_UPDATE;

  // template <>
  // THREADLOCAL mpfr_debug Param<mpfr_debug >::s_epsilon_pivot         = DEFAULT_EPS_PIVOT;

  // Definitions of the Param


  void Param::setEpsilonFactorization(Real eps)
  {
    s_epsilon_factorization = eps;
  }

  Real Param::epsilon()
  {
    return (s_epsilon);
  }

  void Param::setEpsilon(Real eps)
  {
    s_epsilon = eps;
  }


  Real Param::epsilonFactorization()
  {
    return s_epsilon_factorization;
  }

  Real Param::epsilonUpdate()
  {
    return s_epsilon_update;
  }

  void Param::setEpsilonUpdate(Real eps)
  {
    s_epsilon_update = eps;
  }

  Real Param::epsilonPivot()
  {
    return s_epsilon_pivot;
  }

  void Param::setEpsilonPivot(Real eps)
  {
    s_epsilon_pivot = eps;
  }

  bool msginconsistent(const char* name, const char* file, int line)
  {
   assert(name != 0);
   assert(file != 0);
   assert(line >= 0);

   MSG_ERROR( std::cerr << file << "(" << line << ") "
   << "Inconsistency detected in " << name << std::endl; )

   return 0;
  }

  template <>
  Real spxFrexp(Real y, int* exp)
  {
    return frexp(y, exp);
  }

  // @todo: fix this definition
  // template <class T>
  // inline boost::multiprecision::number<T> spxFrexp(boost::multiprecision::number<T> y, int* exp)
  // {
  //   return frexp(y, exp);
  // }

  // @todo: write a boost version of the following function. Check whether this
  // function gets called from the Scalers, if not, we can have a general
  // version of the function in spxdefines.hpp
  Real spxLdexp(Real x, int exp)
  {
    return ldexp(x,exp);
  }

  Rational spxLdexp(Rational x, int exp)
  {
    // This call shouldn't happen. This is a dummy function to deal with the
    // Rational Scalar issue.
    assert(false);
    return 0;
  }




} // namespace soplex
