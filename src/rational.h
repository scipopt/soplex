/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  rational.h
 * @brief Wrapper for GMP types.
 */
#ifndef _RATIONAL_H_
#define _RATIONAL_H_

#include <math.h>
#include <assert.h>
#include <iostream>
#include <string>

#include "spxdefines.h"



#ifdef SOPLEX_WITH_GMP
#include "gmp.h"
#include "gmpxx.h"
#endif

namespace soplex
{
/**@brief   Wrapper for GMP type mpq_class.
 * @ingroup Algebra
 *
 * We wrap mpq_class so that we can replace it by SoPlex's normal Real type if GMP is not available.
 */
#ifdef SOPLEX_WITH_GMP

/// If compiled with GMP support, Rational is defined as mpq_class.
class Rational : public mpq_class
{
public:

   Rational()
      : mpq_class()
   {
   }

   Rational(const Rational& r)
      : mpq_class()
   {
      *this = r;
   }

   Rational(const mpq_class& q)
      : mpq_class(q)
   {
   }

   Rational(const Real& r)
      : mpq_class(r)
   {
   }

   Rational(const int& i)
      : mpq_class(i)
   {
   }

   operator Real() const
   {
      return this->get_d();
   }

   Rational operator+(const Rational& r) const
   {
      Rational retval = *this;
      retval += r;
      return retval;
   }

   Rational operator-(const Rational& r) const
   {
      Rational retval = *this;
      retval -= r;
      return retval;
   }

   Rational operator*(const Rational& r) const
   {
      Rational retval = *this;
      retval *= r;
      return retval;
   }

   Rational operator/(const Rational& r) const
   {
      Rational retval = *this;
      retval /= r;
      return retval;
   }

   /// return as string
   std::string toString(const bool asfloat = true) const;

   /// read Rational from string
   bool readString(const char* s);
};

/// return whether Rational provides exact arithmetic
#define RationalIsExact() (true)

inline static Rational abs(const Rational& r)
{
   Rational res = r;

   if( r < 0 )
      res *= -1;

   return res;
}

/// print Rational
std::ostream& operator<<(std::ostream& os, const Rational& q);

/// Negation.
Rational operator-(const Rational& q);

#else

/// If compiled without GMP support, Rational is defined as SoPlex's normal Real.
typedef Real Rational;

/// return whether Rational provides exact arithmetic
#define RationalIsExact() (false)

#endif

} // namespace soplex
#endif // _RATIONAL_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
