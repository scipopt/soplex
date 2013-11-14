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

namespace soplex
{
/**@brief   Wrapper for GMP type mpq_class.
 * @ingroup Algebra
 *
 * We wrap mpq_class so that we can replace it by a double type if GMP is not available.
 */

/// If compiled with GMP support, Rational is defined as mpq_class.
   class Rational
   {
   private:
      class Private;
      Private* dpointer;

   public:
      /// default constructor
      Rational();

      /// copy constructor
      Rational(const Rational& r);

      /// copy assignment operator
      Rational& operator=(const Rational&);

      /// constructor from long double
      Rational(const long double& r);

      /// constructor from double
      Rational(const double& r);

      ///constructor from int
      Rational(const int& i);

#if 0 /// not currently working

      /// constructor from mpq_class
      Rational(const mpq_class& q);

      /// constructor from mpq_t
      Rational(const mpq_t& q);

#endif

      /// typecasts Rational to double (only allows explicit typecasting)
      explicit operator double() const;

      /// addition operator
      Rational operator+(const Rational& r) const;

      /// addition assignment operator
      Rational operator+=(const Rational& r) const;

      /// addition operator for doubles
      Rational operator+(const double& r) const;

      /// addition assignment operator  for doubles
      Rational operator+=(const double& r) const;

      /// subtraction operator
      Rational operator-(const Rational& r) const;

      /// subtraction assignment operator
      Rational operator-=(const Rational& r) const;

      /// subtraction operator for doubles
      Rational operator-(const double& r) const;

      /// subtraction assignment operator for doubles
      Rational operator-=(const double& r) const;

      /// multiplication operator
      Rational operator*(const Rational& r) const;

      /// multiplication assignment operator operator
      Rational operator*=(const Rational& r) const;

      /// multiplication operator for doubles
      Rational operator*(const double& r) const;

      /// multiplication assignment operator for doubles
      Rational operator*=(const double& r) const;

      /// division operator
      Rational operator/(const Rational& r) const;

      /// division assignment operator
      Rational operator/=(const Rational& r) const;

      /// division operator for doubles
      Rational operator/(const double& r) const;

      /// division assignment operator for doubles
      Rational operator/=(const double& r) const;

      /// read Rational from string
      bool readString(const char* s);

      /// TODO: Place "#define RationalIsExact() (true/false)" in .cpp

      friend std::string rationalToString(const Rational& r, const bool asfloat);
      friend bool readStringRational(const char* s, Rational& value);
      friend Rational abs(const Rational& r);
      friend std::ostream& operator<<(std::ostream& os, const Rational& q);
      friend Rational operator-(const Rational& q);
      friend bool operator!=(const Rational& r, const Rational& s);
      friend bool operator==(const Rational& r, const Rational& s);
      friend bool operator<(const Rational& r, const Rational& s);
      friend bool operator<=(const Rational& r, const Rational& s);
      friend bool operator>(const Rational& r, const Rational& s);
      friend bool operator>=(const Rational& r, const Rational& s);

      friend bool operator!=(const Rational& r, const double& s);
      friend bool operator==(const Rational& r, const double& s);
      friend bool operator<(const Rational& r, const double& s);
      friend bool operator<=(const Rational& r, const double& s);
      friend bool operator>(const Rational& r, const double& s);
      friend bool operator>=(const Rational& r, const double& s);

      friend bool operator!=(const double& r, const Rational& s);
      friend bool operator==(const double& r, const Rational& s);
      friend bool operator<(const double& r, const Rational& s);
      friend bool operator<=(const double& r, const Rational& s);
      friend bool operator>(const double& r, const Rational& s);
      friend bool operator>=(const double& r, const Rational& s);

      friend Rational operator+(const double& d, const Rational& r);
      friend Rational operator-(const double& d, const Rational& r);
      friend Rational operator*(const double& d, const Rational& r);
      friend Rational operator/(const double& d, const Rational& r);
   };

   /// convert rational number to string
   std::string rationalToString(const Rational& r, const bool asfloat = true);

   /// read Rational from string
   bool readStringRational(const char* s, Rational& value);

   /// absolute function
   Rational abs(const Rational& r);

   /// print Rational
   std::ostream& operator<<(std::ostream& os, const Rational& q);

   /// Negation.
   Rational operator-(const Rational& q);

   /// equality operator
   bool operator==(const Rational& r, const Rational& s);

   /// inequality operator
   bool operator!=(const Rational& r, const Rational& s);

   /// less than operator
   bool operator<(const Rational& r, const Rational& s);

   /// less than or equal to operator
   bool operator<=(const Rational& r, const Rational& s);

   /// greater than operator
   bool operator>(const Rational& r, const Rational& s);

   /// greater than or equal to operator
   bool operator>=(const Rational& r, const Rational& s);



   /// equality operator for Rational and double
   bool operator==(const Rational& r, const double& s);

   /// inequality operator for Rational and double
   bool operator!=(const Rational& r, const double& s);

   /// less than operator for Rational and double
   bool operator<(const Rational& r, const double& s);

   /// less than or equal to operator for Rational and double
   bool operator<=(const Rational& r, const double& s);

   /// greater than operator for Rational and double
   bool operator>(const Rational& r, const double& s);

   /// greater than or equal to operator for Rational and double
   bool operator>=(const Rational& r, const double& s);



   /// equality operator for double and Rational
   bool operator==(const double& r, const Rational& s);

   /// inequality operator for double and Rational
   bool operator!=(const double& r, const Rational& s);

   /// less than operator for double and Rational
   bool operator<(const double& r, const Rational& s);

   /// less than or equal to operator for double and Rational
   bool operator<=(const double& r, const Rational& s);

   /// greater than operator for double and Rational
   bool operator>(const double& r, const Rational& s);

   /// greater than or equal to operator for double and Rational
   bool operator>=(const double& r, const Rational& s);


   /// addition operator for double and Rational
   Rational operator+(const double& d, const Rational& r);

   /// addition operator for double and Rational
   Rational operator+(const double& d, const Rational& r);

   /// addition operator for double and Rational
   Rational operator+(const double& d, const Rational& r);

   /// addition operator for double and Rational
   Rational operator+(const double& d, const Rational& r);

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
