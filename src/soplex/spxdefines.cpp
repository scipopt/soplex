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

namespace soplex
{

THREADLOCAL const Real infinity                 = DEFAULT_INFINITY;

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
  template <>
  Real spxLdexp(Real x, int exp)
  {
    return ldexp(x,exp);
  }

  template <>
  Rational spxLdexp(Rational x, int exp)
  {
    // This call shouldn't happen. This is a dummy function to deal with the
    // Rational Scalar issue.
    assert(true);
    return 0;
  }

  template <typename T>
  boost::multiprecision::number<T> spxLdexp(boost::multiprecision::number<T> x, int exp)
  {
    return ldexp(x,exp);
  }



} // namespace soplex
