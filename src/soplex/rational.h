/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2020 Konrad-Zuse-Zentrum                            */
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
#include <string.h>
#include <iostream>

#include "soplex/spxdefines.h"
#include "soplex/idlist.h"


#include "soplex/spxalloc.h"
#ifdef SOPLEX_WITH_BOOST
#include "boost/multiprecision/number.hpp"

#ifdef SOPLEX_WITH_MPFR
#include "boost/multiprecision/mpfr.hpp"
#endif

#ifdef SOPLEX_WITH_CPPMPF
#include "boost/multiprecision/cpp_dec_float.hpp"
#include "boost/multiprecision/cpp_int.hpp"
#endif
#endif

#ifdef SOPLEX_WITH_GMP
#include "gmp.h"
#ifdef SOPLEX_WITH_BOOST
#include "boost/multiprecision/gmp.hpp"
#endif
#endif


namespace soplex
{
/**@brief   Wrapper for GMP type mpq_class.
 * @ingroup Algebra
 *
 * We wrap mpq_class so that we can replace it by a double type if GMP is not available.
 */

/// If compiled with GMP support, Rational is defined as mpq_class.
class Rational // coverity[missing_move_assignment]
{
private:
   class Private;
   Private* dpointer;

#ifdef SOPLEX_WITH_GMP
   THREADLOCAL static IdList< Private > unusedPrivateList;
   THREADLOCAL static bool useListMem;


   /// special constructor only for initializing static rational variables; this is necessary since we need a
   /// constructor for Rational::{ZERO, POSONE, NEGONE} that does not use these numbers
   Rational(const int& i, const bool& dummy);

   ///@name Static variables for special rational values
   ///@{

   static const Rational ZERO;
   static const Rational POSONE;
   static const Rational NEGONE;

   ///@}
#endif

public:

   ///@name Construction and destruction
   ///@{

   /// default constructor
   Rational();

   /// copy constructor
   Rational(const Rational& r);

   /// constructor from long double
   Rational(const long double& r);

   /// constructor from double
   Rational(const double& r);

   ///constructor from int
   Rational(const int& i);

#ifdef SOPLEX_WITH_GMP
   /// constructor from mpq_t (GMP only)
   Rational(const mpq_t& q);
#endif

#ifdef SOPLEX_WITH_BOOST
   // constructor from boost number
   template <typename T, boost::multiprecision::expression_template_option eto>
   Rational(const boost::multiprecision::number<T, eto>& q);
#endif

   /// destructor
   ~Rational();

   /// enables list memory
   static void enableListMem();

   /// frees the unused rational elements in the memory list
   /** this can be useful when you want to save memory or needed when working with a GMP memory manager like the one
    *  in EGlib that frees GMP memory before the destructor of the static memory list is called; in most cases this
    *  method is optional; note that this does not free the Rational elements that are currently in use
    */
   static void freeListMem();

   /// disables list memory
   static void disableListMem();

   /// assignment operator
   Rational& operator=(const Rational&);

   /// assignment operator from long double
   Rational& operator=(const long double& r);

   /// assignment operator from double
   Rational& operator=(const double& r);

   /// assignment operator from int
   Rational& operator=(const int& i);
   ///@}

#ifdef SOPLEX_WITH_GMP
   /// @name GMP Only methods
   ///
   /// Methods of the Rational class that are only available if SoPlex is compiled with "-DGMP=on"
   ///
   ///@{

   /// assignment operator from mpq_t
   Rational& operator=(const mpq_t& q);
#endif

#ifdef SOPLEX_WITH_BOOST
   // assignment operator from boost multiprecision number The operator should
   // convert the boost number to mpq_t

   // Note: the function is only implemented in the #if part
   template <typename T, boost::multiprecision::expression_template_option eto>
   Rational& operator=(const boost::multiprecision::number<T, eto>& q);
   ///@}
#endif

   ///@name Typecasts
   ///@{

   operator double() const;
   operator long double() const;
   operator float() const;
#ifdef SOPLEX_WITH_BOOST
#ifndef SOPLEX_WITH_CPPMPF
   // Operator to typecast Rational to one of the Boost Number types
   template <typename T, boost::multiprecision::expression_template_option eto>
   operator boost::multiprecision::number<T, eto>() const;
#else
   // Operator to typecast Rational to one of the Boost Number types
   template <unsigned bits, boost::multiprecision::expression_template_option eto>
   operator boost::multiprecision::number<boost::multiprecision::backends::cpp_dec_float<bits>, eto>()
   const;
#endif
#endif

#ifdef SOPLEX_WITH_GMP
   /// provides read-only access to underlying mpq_t
   const mpq_t* getMpqPtr() const;

   /// provides read-only access to underlying mpq_t
   const mpq_t& getMpqRef() const;

   /// provides write access to underlying mpq_t; use with care
   mpq_t* getMpqPtr_w() const;

   /// provides write access to underlying mpq_t; use with care
   mpq_t& getMpqRef_w() const;
   ///@} // end of RationalWithGMP
#endif

   ///@name Typecasts
   ///@{

   ///@}


   ///@name Arithmetic operators
   ///@{

   /// addition operator
   Rational operator+(const Rational& r) const;

   /// addition assignment operator
   Rational& operator+=(const Rational& r);

   /// addition operator for doubles
   Rational operator+(const double& r) const;

   /// addition assignment operator  for doubles
   Rational& operator+=(const double& r);

   /// addition operator for ints
   Rational operator+(const int& r) const;

   /// addition assignment operator  for ints
   Rational& operator+=(const int& r);

   /// subtraction operator
   Rational operator-(const Rational& r) const;

   /// subtraction assignment operator
   Rational& operator-=(const Rational& r);

   /// subtraction operator for doubles
   Rational operator-(const double& r) const;

   /// subtraction assignment operator for doubles
   Rational& operator-=(const double& r);

   /// subtraction operator for ints
   Rational operator-(const int& r) const;

   /// subtraction assignment operator for ints
   Rational& operator-=(const int& r);

   /// multiplication operator
   Rational operator*(const Rational& r) const;

   /// multiplication assignment operator operator
   Rational& operator*=(const Rational& r);

   /// multiplication operator for doubles
   Rational operator*(const double& r) const;

   /// multiplication assignment operator for doubles
   Rational& operator*=(const double& r);

   /// multiplication operator for ints
   Rational operator*(const int& r) const;

   /// multiplication assignment operator for ints
   Rational& operator*=(const int& r);

   /// division operator
   Rational operator/(const Rational& r) const;

   /// division assignment operator
   Rational& operator/=(const Rational& r);

   /// division operator for doubles
   Rational operator/(const double& r) const;

   /// division assignment operator for doubles
   Rational& operator/=(const double& r);

   /// division operator for ints
   Rational operator/(const int& r) const;

   /// division assignment operator for ints
   Rational& operator/=(const int& r);

   /// add product of two rationals
   Rational& addProduct(const Rational& r, const Rational& s);

   /// subtract product of two rationals
   Rational& subProduct(const Rational& r, const Rational& s);

   /// add quotient of two rationals, r divided by s
   Rational& addQuotient(const Rational& r, const Rational& s);

   /// subtract quotient of two rationals, r divided by s
   Rational& subQuotient(const Rational& r, const Rational& s);

   /// inversion
   Rational& invert();

   /// round up to next power of two
   Rational& powRound();

   ///@}


   ///@name Methods for checking exactness of doubles
   ///@{

   /// checks if \p d is the closest number that can be represented by double
   bool isNextTo(const double& d);

   /// checks if \p d is exactly equal to the Rational and if not, if it is one of the two adjacent doubles
   bool isAdjacentTo(const double& d) const;

   ///@}


   ///@name Methods for querying size
   ///@{

   /// Size in specified base (bit size for base 2)
   int sizeInBase(const int base = 2) const;

   ///@}


   ///@name Static methods
   ///@{

   /// returns precision of Rational implementation, i.e., number of bits used to store Rational numbers (INT_MAX if exact)
   static int precision();

   ///@}


   ///@name Conversion from and to String
   ///@{

   /// read Rational from string
   bool readString(const char* s);

   friend std::string rationalToString(const Rational& r, const int precision);
   friend bool readStringRational(const char* s, Rational& value);
   friend std::ostream& operator<<(std::ostream& os, const Rational& q);

   ///@}

   ///@name Friends
   ///@{

   friend int compareRational(const Rational& r, const Rational& s);
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

   friend bool operator!=(const Rational& r, const long double& s);
   friend bool operator==(const Rational& r, const long double& s);
   friend bool operator<(const Rational& r, const long double& s);
   friend bool operator<=(const Rational& r, const long double& s);
   friend bool operator>(const Rational& r, const long double& s);
   friend bool operator>=(const Rational& r, const long double& s);

   friend bool operator!=(const long double& r, const Rational& s);
   friend bool operator==(const long double& r, const Rational& s);
   friend bool operator<(const long double& r, const Rational& s);
   friend bool operator<=(const long double& r, const Rational& s);
   friend bool operator>(const long double& r, const Rational& s);
   friend bool operator>=(const long double& r, const Rational& s);

   friend bool operator!=(const Rational& r, const float& s);
   friend bool operator==(const Rational& r, const float& s);
   friend bool operator<(const Rational& r, const float& s);
   friend bool operator<=(const Rational& r, const float& s);
   friend bool operator>(const Rational& r, const float& s);
   friend bool operator>=(const Rational& r, const float& s);

   friend bool operator!=(const float& r, const Rational& s);
   friend bool operator==(const float& r, const Rational& s);
   friend bool operator<(const float& r, const Rational& s);
   friend bool operator<=(const float& r, const Rational& s);
   friend bool operator>(const float& r, const Rational& s);
   friend bool operator>=(const float& r, const Rational& s);


   friend Rational operator+(const double& d, const Rational& r);
   friend Rational operator-(const double& d, const Rational& r);
   friend Rational operator*(const double& d, const Rational& r);
   friend Rational operator/(const double& d, const Rational& r);

   friend bool operator!=(const Rational& r, const int& s);
   friend bool operator==(const Rational& r, const int& s);
   friend bool operator<(const Rational& r, const int& s);
   friend bool operator<=(const Rational& r, const int& s);
   friend bool operator>(const Rational& r, const int& s);
   friend bool operator>=(const Rational& r, const int& s);

   friend bool operator!=(const int& r, const Rational& s);
   friend bool operator==(const int& r, const Rational& s);
   friend bool operator<(const int& r, const Rational& s);
   friend bool operator<=(const int& r, const Rational& s);
   friend bool operator>(const int& r, const Rational& s);
   friend bool operator>=(const int& r, const Rational& s);

   friend Rational operator+(const int& d, const Rational& r);
   friend Rational operator-(const int& d, const Rational& r);
   friend Rational operator*(const int& d, const Rational& r);
   friend Rational operator/(const int& d, const Rational& r);

   friend Rational spxAbs(const Rational& r);
   friend int sign(const Rational& r);
   friend Rational operator-(const Rational& q);

   ///@}
};

/// less than operator
bool operator<(const Rational& r, const Rational& s);

///@name Parsing and printing
///@{

/// convert rational number to string
std::string rationalToString(const Rational& r, const int precision = 32);

/// read Rational from string
bool readStringRational(const char* s, Rational& value);

/// print Rational
std::ostream& operator<<(std::ostream& os, const Rational& r);

///@}

/// less than operator for Rational and double
bool operator<(const Rational& r, const double& s);

///@name Relational operators
///@{

/// comparison operator returning a positive value if r > s, zero if r = s, and a negative value if r < s
int compareRational(const Rational& r, const Rational& s);

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

/// equality operator for Rational and long double
bool operator==(const Rational& r, const long double& s);

/// inequality operator for Rational and long double
bool operator!=(const Rational& r, const long double& s);

/// less than operator for Rational and long double
bool operator<(const Rational& r, const long double& s);

/// less than or equal to operator for Rational and long double
bool operator<=(const Rational& r, const long double& s);

/// greater than operator for Rational and long double
bool operator>(const Rational& r, const long double& s);

/// greater than or equal to operator for Rational and long double
bool operator>=(const Rational& r, const long double& s);

/// equality operator for long double and Rational
bool operator==(const long double& r, const Rational& s);

/// inequality operator for long double and Rational
bool operator!=(const long double& r, const Rational& s);

/// less than operator for long double and Rational
bool operator<(const long double& r, const Rational& s);

/// less than or equal to operator for long double and Rational
bool operator<=(const long double& r, const Rational& s);

/// greater than operator for long double and Rational
bool operator>(const long double& r, const Rational& s);

/// greater than or equal to operator for long double and Rational
bool operator>=(const long double& r, const Rational& s);

/// equality operator for Rational and int
bool operator==(const Rational& r, const int& s);

/// inequality operator for Rational and int
bool operator!=(const Rational& r, const int& s);

/// less than operator for Rational and int
bool operator<(const Rational& r, const int& s);

/// less than or equal to operator for Rational and int
bool operator<=(const Rational& r, const int& s);

/// greater than operator for Rational and int
bool operator>(const Rational& r, const int& s);

/// greater than or equal to operator for Rational and int
bool operator>=(const Rational& r, const int& s);

/// equality operator for int and Rational
bool operator==(const int& r, const Rational& s);

/// inequality operator for int and Rational
bool operator!=(const int& r, const Rational& s);

/// less than operator for int and Rational
bool operator<(const int& r, const Rational& s);

/// less than or equal to operator for int and Rational
bool operator<=(const int& r, const Rational& s);

/// greater than operator for int and Rational
bool operator>(const int& r, const Rational& s);

/// greater than or equal to operator for int and Rational
bool operator>=(const int& r, const Rational& s);

///@}

/// Sign function; returns 1 if r > 0, 0 if r = 0, and -1 if r < 0.
int sign(const Rational& r);

///@name Non-member arithmetic operators and functions
///@{

/// addition operator for double and Rational
Rational operator+(const double& d, const Rational& r);

/// addition operator for double and Rational
Rational operator+(const double& d, const Rational& r);

/// addition operator for double and Rational
Rational operator+(const double& d, const Rational& r);

/// addition operator for double and Rational
Rational operator+(const double& d, const Rational& r);

/// addition operator for int and Rational
Rational operator+(const int& d, const Rational& r);

/// addition operator for int and Rational
Rational operator+(const int& d, const Rational& r);

/// addition operator for int and Rational
Rational operator+(const int& d, const Rational& r);

/// addition operator for int and Rational
Rational operator+(const int& d, const Rational& r);

/// absolute function
Rational spxAbs(const Rational& r);

/// Sign function; returns 1 if r > 0, 0 if r = 0, and -1 if r < 0.
int sign(const Rational& r);

/// Negation.
Rational operator-(const Rational& r);

/// Total size of rational vector.
int totalSizeRational(const Rational* vector, const int length, const int base = 2);

/// Size of least common multiple of denominators in rational vector.
int dlcmSizeRational(const Rational* vector, const int length, const int base = 2);

/// Size of largest denominator in rational vector.
int dmaxSizeRational(const Rational* vector, const int length, const int base = 2);

#ifdef SOPLEX_WITH_GMP
/// Defines the "Pimpl"-class Private
class Rational::Private
{
public:

   mpq_t privatevalue;  ///< actual value of the Rational object
   Private* theprev;    ///< pointer to the previous element in the list
   Private* thenext;    ///< pointer to the next element in the list

   /// default constructor
   Private()
      : theprev(0)
      , thenext(0)
   {
      mpq_init(privatevalue);
   }

   /// copy constructor
   Private(const Private& p)
      : theprev(0)
      , thenext(0)
   {
      // a newly constructed element is not in any list, even if the original element (p) is; hence we initialize
      // theprev and thenext to zero
      mpq_init(privatevalue);
      mpq_set(this->privatevalue, p.privatevalue);
   }

   /// constructor from long double
   Private(const long double& r)
      : theprev(0)
      , thenext(0)
   {
      mpq_init(privatevalue);

      if(r == (long double)(1.0))
         mpq_set(privatevalue, Rational::POSONE.dpointer->privatevalue);
      else if(r == (long double)(-1.0))
         mpq_set(privatevalue, Rational::NEGONE.dpointer->privatevalue);
      else if(r == (long double)(0.0))
      {
         assert(mpq_equal(privatevalue, Rational::ZERO.dpointer->privatevalue) != 0);
      }
      else
         mpq_set_d(privatevalue, double(r));
   }

   /// constructor from double
   Private(const double& r)
      : theprev(0)
      , thenext(0)
   {
      mpq_init(privatevalue);

      if(r == 1.0)
         mpq_set(privatevalue, Rational::POSONE.dpointer->privatevalue);
      else if(r == -1.0)
         mpq_set(privatevalue, Rational::NEGONE.dpointer->privatevalue);
      else if(r == 0.0)
      {
         assert(mpq_equal(privatevalue, Rational::ZERO.dpointer->privatevalue) != 0);
      }
      else
         mpq_set_d(privatevalue, r);
   }

   /// constructor from int
   Private(const int& i)
      : theprev(0)
      , thenext(0)
   {
      mpq_init(privatevalue);

      if(i == 1)
         mpq_set(privatevalue, Rational::POSONE.dpointer->privatevalue);
      else if(i == -1)
         mpq_set(privatevalue, Rational::NEGONE.dpointer->privatevalue);
      else if(i == 0)
      {
         assert(mpq_equal(privatevalue, Rational::ZERO.dpointer->privatevalue) != 0);
      }
      else
         mpq_set_si(privatevalue, i, 1);
   }

   /// constructor from mpq_t (GMP only)
   Private(const mpq_t& q)
      : theprev(0)
      , thenext(0)
   {
      mpq_init(privatevalue);
      mpq_set(privatevalue, q);
   }

   /// constructor from boost number
   // Also should satisfy SOPLEX_WITH_GMP
#ifdef SOPLEX_WITH_BOOST
#ifdef SOPLEX_WITH_GMP
   template <typename T, boost::multiprecision::expression_template_option eto>
   Private(const boost::multiprecision::number<T, eto>& q)
      : theprev(0)
      , thenext(0)
   {
      boost::multiprecision::mpq_rational tmp{q};
      mpq_init(privatevalue);
      mpq_set(privatevalue, tmp.backend().data());
   }
#endif
#endif

   /// destructor
   ~Private()
   {
      mpq_clear(privatevalue);
   }

   /// assignment operator
   Private& operator=(const Private& p)
   {
#ifdef SOPLEX_PERFALT_4

      if(mpq_equal(this->privatevalue, p.privatevalue) != 0)
         return *this;

#endif

      // we only assign the value; the position in the list, i.e., theprev and thenext, must not be modified
      mpq_set(this->privatevalue, p.privatevalue);
      return *this;
   }

   /// assignment operator from long double
   Private& operator=(const long double& r)
   {
      // we only assign the value; the position in the list, i.e., theprev and thenext, must not be modified
      if(r == (long double)(0.0))
      {
#ifdef SOPLEX_PERFALT_5a
#ifdef SOPLEX_PERFALT_1

         if(mpq_sgn(privatevalue) != 0)
#else
         if(mpq_equal(privatevalue, Rational::ZERO.dpointer->privatevalue) == 0)
#endif
#endif
            mpq_set(privatevalue, Rational::ZERO.dpointer->privatevalue);
      }
      else if(r == (long double)(1.0))
      {
#ifdef SOPLEX_PERFALT_5b

         if(mpq_equal(privatevalue, Rational::POSONE.dpointer->privatevalue) == 0)
#endif
            mpq_set(privatevalue, Rational::POSONE.dpointer->privatevalue);
      }
      else if(r == (long double)(-1.0))
      {
#ifdef SOPLEX_PERFALT_5b

         if(mpq_equal(privatevalue, Rational::NEGONE.dpointer->privatevalue) == 0)
#endif
            mpq_set(privatevalue, Rational::NEGONE.dpointer->privatevalue);
      }
      else
         mpq_set_d(this->privatevalue, double(r));

      return *this;
   }

   /// assignment operator from double
   Private& operator=(const double& r)
   {
      // we only assign the value; the position in the list, i.e., theprev and thenext, must not be modified
      if(r == 0.0)
      {
#ifdef SOPLEX_PERFALT_5a
#ifdef SOPLEX_PERFALT_1

         if(mpq_sgn(privatevalue) != 0)
#else
         if(mpq_equal(privatevalue, Rational::ZERO.dpointer->privatevalue) == 0)
#endif
#endif
            mpq_set(privatevalue, Rational::ZERO.dpointer->privatevalue);
      }
      else if(r == 1.0)
      {
#ifdef SOPLEX_PERFALT_5b

         if(mpq_equal(privatevalue, Rational::POSONE.dpointer->privatevalue) == 0)
#endif
            mpq_set(privatevalue, Rational::POSONE.dpointer->privatevalue);
      }
      else if(r == -1.0)
      {
#ifdef SOPLEX_PERFALT_5b

         if(mpq_equal(privatevalue, Rational::NEGONE.dpointer->privatevalue) == 0)
#endif
            mpq_set(privatevalue, Rational::NEGONE.dpointer->privatevalue);
      }
      else
         mpq_set_d(privatevalue, r);

      return *this;
   }

   /// assignment operator from int
   Private& operator=(const int& i)
   {
      // we only assign the value; the position in the list, i.e., theprev and thenext, must not be modified
      if(i == 0)
      {
#ifdef SOPLEX_PERFALT_5a
#ifdef SOPLEX_PERFALT_1

         if(mpq_sgn(privatevalue) != 0)
#else
         if(mpq_equal(privatevalue, Rational::ZERO.dpointer->privatevalue) == 0)
#endif
#endif
            mpq_set(privatevalue, Rational::ZERO.dpointer->privatevalue);
      }
      else if(i == 1)
      {
#ifdef SOPLEX_PERFALT_5b

         if(mpq_equal(privatevalue, Rational::POSONE.dpointer->privatevalue) == 0)
#endif
            mpq_set(privatevalue, Rational::POSONE.dpointer->privatevalue);
      }
      else if(i == -1)
      {
#ifdef SOPLEX_PERFALT_5b

         if(mpq_equal(privatevalue, Rational::NEGONE.dpointer->privatevalue) == 0)
#endif
            mpq_set(privatevalue, Rational::NEGONE.dpointer->privatevalue);
      }
      else
         mpq_set_si(privatevalue, i, 1);

      return *this;
   }

   /// assignment operator from mpq_t
   Private& operator=(const mpq_t& q)
   {
#ifdef SOPLEX_PERFALT_4

      if(mpq_equal(this->privatevalue, q) != 0)
         return *this;

#endif

      // we only assign the value; the position in the list, i.e., theprev and thenext, must not be modified
      mpq_set(this->privatevalue, q);
      return *this;
   }

#ifdef SOPLEX_WITH_BOOST
   // The back end for Rational operator=.
   template <typename T, boost::multiprecision::expression_template_option eto>
   Private& operator=(const boost::multiprecision::number<T, eto>& q)
   {
      // mpq_rational is used as an intermediate to convert the mpf float to
      // mpq_t.
      boost::multiprecision::mpq_rational tmp{q};
      mpq_set(this->privatevalue, tmp.backend().data());
      return *this;
   }
#endif

   /// previous Private element
   Private*& prev()
   {
      return theprev;
   }

   /// previous Private element
   Private* const& prev() const
   {
      return theprev;
   }

   /// next Private element
   Private*& next()
   {
      return thenext;
   }

   /// next Private element
   Private* const& next() const
   {
      return thenext;
   }
};
#else  // SOPLEX_WITH_GMP

/// Defines the "Pimpl"-class Private
class Rational::Private
{

public:

   /// value
   long double privatevalue;

#ifdef SOPLEX_WITH_BOOST
   template <typename T, boost::multiprecision::expression_template_option eto>
   Private(const boost::multiprecision::number<T, eto>& q)
   {
      privatevalue = (long double)q;
   }

   template <typename T, boost::multiprecision::expression_template_option eto>
   Private& operator=(const boost::multiprecision::number<T, eto>& q)
   {
      privatevalue = (long double)q;
      return *this;
   }
#endif


   /// default constructor
   Private()
   {
      privatevalue = 0;
   }

   /// copy constructor
   Private(const Private& p)
   {
      *this = p;
   }

   /// constructor from long double
   Private(const long double& r)
   {
      privatevalue = r;
   }

   /// constructor from double
   Private(const double& r)
   {
      privatevalue = r;
   }

   /// constructor from int
   Private(const int& i)
   {
      privatevalue = i;
   }

   /// assignment operator
   Private& operator=(const Private& p)
   {
      this->privatevalue = p.privatevalue;
      return *this;
   }

   /// assignment operator from long double
   Private& operator=(const long double& r)
   {
      this->privatevalue = r;
      return *this;
   }

   /// assignment operator from double
   Private& operator=(const double& r)
   {
      this->privatevalue = (long double)(r);
      return *this;
   }

   /// assignment operator from int
   Private& operator=(const int& i)
   {
      this->privatevalue = (long double)(i);
      return *this;
   }
};
#endif // SOPLEX_WITH_GMP



#ifdef SOPLEX_WITH_BOOST
using namespace boost::multiprecision;
// Definitions related to boost number and SoPlex Rational. This is still
// inside the #ifdef SOPLEX_WITH_GMP

// Assignment operator from boost number. Uses the API for the Private class
// to do this.
template <typename T, expression_template_option eto>
Rational& Rational::operator=(const number<T, eto>& r)
{
   *(this->dpointer) = r;
   return *this;
}

#ifndef SOPLEX_WITH_CPPMPF
// Operator to typecast Rational to one of the Boost Number types
template <typename T, expression_template_option eto>
Rational::operator number<T,  eto>() const
{
   // Constructs a boost::multiprecision::number<T> with value
   // this->pointer->privatevalue
   return number<T, eto>(this->dpointer->privatevalue);
}

#else
#ifdef SOPLEX_WITH_GMP
// Specialization for the conversion mpq_t -> cpp_rational
template<unsigned bits, expression_template_option eto>
Rational::operator number<backends::cpp_dec_float<bits>, eto>() const
{
   number<gmp_rational, et_on> mpq_numb(this->dpointer->privatevalue);
   number<cpp_rational_backend, et_on> cpp_numb = cpp_rational(mpq_numb);
   return number<backends::cpp_dec_float<bits>, eto>(cpp_numb);
}
#else
// Specialization for the conversion double -> cpp_rational
template<unsigned bits, expression_template_option eto>
Rational::operator number<backends::cpp_dec_float<bits>, eto>() const
{
   return number<backends::cpp_dec_float<bits>, eto>(this->dpointer->privatevalue);
}
#endif
#endif

// Constructor from boost number. Code is exactly same as that of construction
// of Rational from mpq_t, basically a wrapper around the assignment operator
// (=) of the Private class; calls the boost number assignment operator.
#ifdef SOPLEX_WITH_MPFR
template <typename T, boost::multiprecision::expression_template_option eto>
Rational::Rational(const boost::multiprecision::number<T, eto>& q)
{

   // TODO Figure out why SCIP complains about the static variable problem
   // (that useListMem doesn't exist)
   if(Rational::useListMem)
   {
      dpointer = unusedPrivateList.last();

      if(dpointer != nullptr)
      {
         assert(unusedPrivateList.first() != 0);
         unusedPrivateList.remove(dpointer);
         *dpointer = q;
      }
      else
      {
         assert(unusedPrivateList.first() == 0);
         spx_alloc(dpointer);
         new(dpointer) Private(q);
      }
   }
   else
   {
      assert(unusedPrivateList.length() == 0);
      dpointer = 0;
      spx_alloc(dpointer);
      new(dpointer) Private(q);
   }

   assert(dpointer != 0);

}
#endif  // SOPLEX_WITH_CPPMPF

#ifdef SOPLEX_WITH_CPPMPF
template <typename T, boost::multiprecision::expression_template_option eto>
Rational::Rational(const boost::multiprecision::number<T, eto>& q)
{
   dpointer = 0;
   spx_alloc(dpointer);
   new(dpointer) Private(q);

   assert(dpointer != 0);
}
#endif // SOPLEX_WITH_CPPMPF
#endif

// A dummy function to deal with the rational scalar issue. This will never be
// called
inline Rational spxFrexp(Rational r, int* d)
{
   assert(false);
   return Rational(0);
}

// same as before
inline Rational spxLdexp(Rational x, int exp)
{
   // This call shouldn't happen. This is a dummy function to deal with the
   // Rational Scalar issue.
   assert(false);
   return 0;
}



//@}

} // namespace soplex

#endif // _RATIONAL_H_
