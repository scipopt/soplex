/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  rational.cpp
 * @brief Wrapper for GMP types.
 */

#ifndef SOPLEX_LEGACY
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <limits.h>


#include "rational.h"
#include "spxalloc.h"
#include "spxdefines.h"

#ifdef SOPLEX_WITH_GMP
#include "gmp.h"
#endif

#ifdef SOPLEX_WITH_GMPXX
#include "gmp.h"
#include "gmpxx.h"
#endif

namespace soplex
{
#ifdef SOPLEX_WITH_GMP

/// Defines the "Pimpl"-class Private
class Rational::Private
{
public:

   /// actual value of the Rational object
   mpq_t privatevalue;

   /// default constructor
   Private()
   {
      mpq_init(privatevalue);
   }

   /// copy constructor
   Private(const Private& p)
   {
      mpq_init(privatevalue);
      mpq_set(this->privatevalue, p.privatevalue);
   }

   /// copy assignment operator
   Private& operator=(const Private& p)
   {
      mpq_set(this->privatevalue, p.privatevalue);
      return *this;
   }

   /// constructor from long double
   Private(const long double& r)
   {
      mpq_init(privatevalue);
      mpq_set_d(privatevalue, double(r));
   }

   /// constructor from double
   Private(const double& r)
   {
      mpq_init(privatevalue);
      mpq_set_d(privatevalue, r);
   }

   /// constructor from int
   Private(const int& i)
   {
      mpq_init(privatevalue);
      mpq_set_d(privatevalue, i);
   }

   /// constructor from mpq_t
   Private(const mpq_t& q)
   {
      mpq_init(privatevalue);
      mpq_set(privatevalue, q);
   }

   /// destructor
   ~Private()
   {
      mpq_clear(privatevalue);
   }
};



/// default constructor
Rational::Rational()
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private();
}



/// copy constructor
Rational::Rational(const Rational& r)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(*(r.dpointer));
}



/// copy assignment operator
Rational& Rational::operator=(const Rational &r)
{
   *(this->dpointer) = *(r.dpointer);
   return *this;
}



/// constructor from long double
Rational::Rational(const long double& r)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(r);
}



/// constructor from double
Rational::Rational(const double& r)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(r);
}



/// constructor from int
Rational::Rational(const int& i)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(i);
}



/// constructor from mpq_t
Rational::Rational(const mpq_t& q)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(q);
}



/// destructor
Rational::~Rational()
{
   dpointer->~Private();
   spx_free(dpointer);
}



/// typecasts Rational to double (allows only explicit typecast)
Rational::operator double() const
{
   return mpq_get_d(this->dpointer->privatevalue);
}



/// typecasts Rational to long double (allows only explicit typecast)
Rational::operator long double() const
{
   return (long double)mpq_get_d(this->dpointer->privatevalue);
}



/// addition operator
Rational Rational::operator+(const Rational& r) const
{
   Rational retval;
   mpq_add(retval.dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return retval;
}



/// addition assignment operator
Rational Rational::operator+=(const Rational& r) const
{
   mpq_add(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return *this;
}



/// addition operator for doubles
Rational Rational::operator+(const double& d) const
{
   Rational retval(d);
   mpq_add(retval.dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
   return retval;
}



/// addition assignment operator for doubles
Rational Rational::operator+=(const double& d) const
{
   mpq_t doubleval;
   mpq_init(doubleval);
   mpq_set_d(doubleval, d);
   mpq_add(this->dpointer->privatevalue, this->dpointer->privatevalue, doubleval);
   mpq_clear(doubleval);
   return *this;
}



/// subtraction operator
Rational Rational::operator-(const Rational& r) const
{
   Rational retval;
   mpq_sub(retval.dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return retval;
}



/// subtraction assignment operator
Rational Rational::operator-=(const Rational& r) const
{
   mpq_sub(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return *this;
}



/// subtraction operator for doubles
Rational Rational::operator-(const double& d) const
{
   Rational retval(d);
   mpq_sub(retval.dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
   return retval;
}



/// subtraction assignment operator for doubles
Rational Rational::operator-=(const double& d) const
{
   mpq_t doubleval;
   mpq_init(doubleval);
   mpq_set_d(doubleval, d);
   mpq_sub(this->dpointer->privatevalue, this->dpointer->privatevalue, doubleval);
   mpq_clear(doubleval);
   return *this;
}



/// multiplication operator
Rational Rational::operator*(const Rational& r) const
{
   Rational retval;
   mpq_mul(retval.dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return retval;
}



/// multiplication assignment operator
Rational Rational::operator*=(const Rational& r) const
{
   mpq_mul(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return *this;
}



/// multiplication operator for doubles
Rational Rational::operator*(const double& d) const
{
   Rational retval(d);
   mpq_mul(retval.dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
   return retval;
}



/// multiplication assignment operator for doubles
Rational Rational::operator*=(const double& d) const
{
   mpq_t doubleval;
   mpq_init(doubleval);
   mpq_set_d(doubleval, d);
   mpq_mul(this->dpointer->privatevalue, this->dpointer->privatevalue, doubleval);
   mpq_clear(doubleval);
   return *this;
}



/// division operator
Rational Rational::operator/(const Rational& r) const
{
   Rational retval;
   mpq_div(retval.dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return retval;
}



/// division assignment operator
Rational Rational::operator/=(const Rational& r) const
{
   mpq_div(this->dpointer->privatevalue, this->dpointer->privatevalue, r.dpointer->privatevalue);
   return *this;
}



/// division operator for doubles
Rational Rational::operator/(const double& d) const
{
   Rational retval(d);
   mpq_div(retval.dpointer->privatevalue, this->dpointer->privatevalue, retval.dpointer->privatevalue);
   return retval;
}



/// division assignment operator for doubles
Rational Rational::operator/=(const double& d) const
{
   mpq_t doubleval;
   mpq_init(doubleval);
   mpq_set_d(doubleval, d);
   mpq_div(this->dpointer->privatevalue, this->dpointer->privatevalue, doubleval);
   mpq_clear(doubleval);
   return *this;
}



/// checks if d is the closest possible double
bool Rational::isNextTo(const double& d)
{
   // get intervall [a,b] of doubles that the Rational is in
   double x = mpq_get_d(this->dpointer->privatevalue);
   double a;
   double b;

   if( Rational(x) < *this )
   {
      a = x;
      b = nextafter(a, infinity);
   }
   else
   {
      b = x;
      a = nextafter(b, -infinity);
   }

   // check if d equals the closer end of the intervall
   bool result = (abs(*this - a) < abs(*this - b))
      ? (d == a)
      : (d == b);

   return result;
}



/// checks if d is exactly equal to the Rational and if not, if it is one of the two adjacent doubles
bool Rational::isAdjacentTo(const double& d)
{
   double x = mpq_get_d(this->dpointer->privatevalue);
   double a;
   double b;
   mpq_t tmp;

   mpq_init(tmp);
   mpq_set_d(tmp, x);

   int cmp = mpq_cmp(tmp, this->dpointer->privatevalue);
   mpq_clear(tmp);

   // the rounded value is smaller than the rational value
   if( cmp < 0 )
   {
      a = x;
      b = nextafter(a, infinity);
   }
   // the rounded value is larger than the rational value
   else if( cmp > 0 )
   {
      b = x;
      a = nextafter(b, -infinity);
   }
   // the rational value is representable in double precision
   else
      return (x == d);

   return ((a == d) || (b == d));
}



/// returns precision of Rational implementation, i.e., number of bits used to store Rational numbers (INT_MAX if exact)
int Rational::precision()
{
   return INT_MAX;
}



#define MAX_STR_LEN 10000
/// read Rational from string
bool Rational::readString(const char* s)
{
   assert(s != 0);
   assert(strlen(s) <= MAX_STR_LEN);
   Rational value;
   const char* pos;

   // if there is a slash or there is no dot and exponent (i.e. we
   // have an integer), we may simply call GMP's string reader
   if( strchr(s, '/') != 0 || strpbrk(s, ".eE") == 0 )
   {
      pos = (*s == '+') ? s + 1 : s;
      if( mpq_set_str(value.dpointer->privatevalue, pos, 10) == 0 )
      {
         mpq_canonicalize(value.dpointer->privatevalue);
         mpq_set(this->dpointer->privatevalue, value.dpointer->privatevalue);
         return true;
      }
      else
         return false;
   }

   // otherwise we analyze the string
#ifndef NDEBUG
   bool has_exponent = false;
   bool has_dot = false;
#endif
   bool has_digits = false;
   bool has_emptyexponent = false;
   long int exponent = 0;
   long int decshift = 0;
   mpz_t shiftpower;
   mpz_init(shiftpower);
   mpq_t shiftpowerRational;
   mpq_init(shiftpowerRational);
   char* t;
   char tmp[MAX_STR_LEN];

   pos = s;

   // 1. sign
   if( (*pos == '+') || (*pos == '-') )
      pos++;

   // 2. Digits before the decimal dot
   while( (*pos >= '0') && (*pos <= '9') )
   {
      has_digits = true;
      pos++;
   }

   // 3. Decimal dot
   if( *pos == '.' )
   {
#ifndef NDEBUG
      has_dot = true;
#endif
      pos++;

      // 4. If there was a dot, possible digit behind it
      while( (*pos >= '0') && (*pos <= '9') )
      {
         has_digits = true;
         pos++;
      }
   }

   // 5. Exponent
   if( tolower(*pos) == 'e' )
   {
#ifndef NDEBUG
      has_exponent = true;
#endif
      has_emptyexponent = true;
      pos++;

      // 6. Exponent sign
      if( (*pos == '+') || (*pos == '-') )
         pos++;

      // 7. Exponent digits
      while( (*pos >= '0') && (*pos <= '9') )
      {
         has_emptyexponent = false;
         pos++;
      }
   }

   if( has_emptyexponent || !has_digits )
      return false;

   assert(has_exponent || has_dot);

   //read up to dot recording digits
   t = tmp;
   pos = s;

   if( *pos == '+' )
      pos++;

   while( ((*pos >= '0') && (*pos <= '9') ) || *pos == '+' || *pos == '-'  )
   {
      *t++ = *pos;
      pos++;
   }
   //record digits after dot, recording positions
   decshift = 0;
   if( *pos == '.' )
   {
      assert(has_dot);
      pos++;
      while( (*pos >= '0') && (*pos <= '9') )
      {
         *t++ = *pos;
         decshift++;
         pos++;
      }
   }
   *t = '\0';

   if( mpq_set_str(value.dpointer->privatevalue, tmp, 10) != 0)
      return false;
   mpq_canonicalize(value.dpointer->privatevalue);

   //record exponent and update final result
   exponent = -decshift;
   if( tolower(*pos) == 'e' )
   {
      pos++;
      assert(has_exponent);
      for( t = tmp; *pos != '\0'; pos++ )
         *t++ = *pos;
      *t = '\0';
      exponent += atol(tmp);
   }
   if( exponent > 0 )
   {
      mpz_ui_pow_ui(shiftpower, 10, exponent);
      mpq_set_z(shiftpowerRational, shiftpower);
      mpq_mul(value.dpointer->privatevalue, value.dpointer->privatevalue, shiftpowerRational);
   }
   else if( exponent < 0 )
   {
      mpz_ui_pow_ui(shiftpower, 10, -exponent);
      mpq_set_z(shiftpowerRational, shiftpower);
      mpq_div(value.dpointer->privatevalue, value.dpointer->privatevalue, shiftpowerRational);
   }

   mpq_canonicalize(value.dpointer->privatevalue);
   mpq_set(this->dpointer->privatevalue, value.dpointer->privatevalue);
   mpz_clear(shiftpower);
   mpq_clear(shiftpowerRational);
   return true;

}



/// convert rational number to string
std::string rationalToString(const Rational& r, const bool asfloat)
{
   std::stringstream sstream;

   sstream << r;
   if( !asfloat )
      return sstream.str();
   else
   {
      sstream.str("");
      sstream << std::setprecision(20) << double(r);
      return sstream.str();
   }

   ///@todo the following code creates a string of the form 1.333...e-06; this would be nice for human readable output,
   ///      because it indicates that precision may have been lost, however, this cannot be parsed by our awk evaluation
   ///      scripts; we should think about how this could be useful
#if 0
   sstream.str("");
   sstream << mpf_class(r);
   std::string origstring = sstream.str();
   sstream.str("");

   size_t epos = origstring.find("e");

   sstream << origstring.substr(0, epos);

   if( origstring.find(".") == std::string::npos )
      sstream << ".0";

   sstream << "...";

   if( epos != std::string::npos )
      sstream << origstring.substr(origstring.find("e"));

   return sstream.str();
#endif
}



/// read Rational from string
bool readStringRational(const char* s, Rational& value)
{
   assert(s != 0);
   assert(strlen(s) <= MAX_STR_LEN);
   const char* pos;

   // if there is a slash or there is no dot and exponent (i.e. we
   // have an integer), we may simply call GMP's string reader
   if( strchr(s, '/') != 0 || strpbrk(s, ".eE") == 0 )
   {
      pos = (*s == '+') ? s + 1 : s;
      if( mpq_set_str(value.dpointer->privatevalue, pos, 10) == 0 )
      {
         mpq_canonicalize(value.dpointer->privatevalue);
         return true;
      }
      else
         return false;
   }

   // otherwise we analyze the string
#ifndef NDEBUG
   bool has_exponent = false;
   bool has_dot = false;
#endif
   bool has_digits = false;
   bool has_emptyexponent = false;
   long int exponent = 0;
   long int decshift = 0;
   mpz_t shiftpower;
   mpz_init(shiftpower);
   mpq_t shiftpowerRational;
   mpq_init(shiftpowerRational);
   char* t;
   char tmp[MAX_STR_LEN];

   pos = s;

   // 1. sign
   if( (*pos == '+') || (*pos == '-') )
      pos++;

   // 2. Digits before the decimal dot
   while( (*pos >= '0') && (*pos <= '9') )
   {
      has_digits = true;
      pos++;
   }

   // 3. Decimal dot
   if( *pos == '.' )
   {
#ifndef NDEBUG
      has_dot = true;
#endif
      pos++;

      // 4. If there was a dot, possible digit behind it
      while( (*pos >= '0') && (*pos <= '9') )
      {
         has_digits = true;
         pos++;
      }
   }

   // 5. Exponent
   if( tolower(*pos) == 'e' )
   {
#ifndef NDEBUG
      has_exponent = true;
#endif
      has_emptyexponent = true;
      pos++;

      // 6. Exponent sign
      if( (*pos == '+') || (*pos == '-') )
         pos++;

      // 7. Exponent digits
      while( (*pos >= '0') && (*pos <= '9') )
      {
         has_emptyexponent = false;
         pos++;
      }
   }

   if( has_emptyexponent || !has_digits )
      return false;

   assert(has_exponent || has_dot);

   // read up to dot recording digits
   t = tmp;
   pos = s;

   if( *pos == '+' )
      pos++;

   while( ((*pos >= '0') && (*pos <= '9') ) || *pos == '+' || *pos == '-'  )
   {
      *t++ = *pos;
      pos++;
   }
   // record digits after dot, recording positions
   decshift = 0;
   if( *pos == '.' )
   {
      assert(has_dot);
      pos++;
      while( (*pos >= '0') && (*pos <= '9') )
      {
         *t++ = *pos;
         decshift++;
         pos++;
      }
   }
   *t = '\0';

   if( mpq_set_str(value.dpointer->privatevalue, tmp, 10) != 0)
      return false;

   mpq_canonicalize(value.dpointer->privatevalue);

   // record exponent and update final result
   exponent = -decshift;
   if( tolower(*pos) == 'e' )
   {
      pos++;
      assert(has_exponent);
      for( t = tmp; *pos != '\0'; pos++ )
         *t++ = *pos;
      *t = '\0';
      exponent += atol(tmp);
   }
   if( exponent > 0 )
   {
      mpz_ui_pow_ui(shiftpower, 10, exponent);
      mpq_set_z(shiftpowerRational, shiftpower);
      mpq_mul(value.dpointer->privatevalue, value.dpointer->privatevalue, shiftpowerRational);
   }
   else if( exponent < 0 )
   {
      mpz_ui_pow_ui(shiftpower, 10, -exponent);
      mpq_set_z(shiftpowerRational, shiftpower);
      mpq_div(value.dpointer->privatevalue, value.dpointer->privatevalue, shiftpowerRational);
   }

   mpq_canonicalize(value.dpointer->privatevalue);
   mpz_clear(shiftpower);
   mpq_clear(shiftpowerRational);

   return true;
}



/// print Rational
std::ostream& operator<<(std::ostream& os, const Rational& r)
{
   char* buffer;
   buffer = (char*) malloc (mpz_sizeinbase(mpq_numref(r.dpointer->privatevalue), 10) + mpz_sizeinbase(mpq_denref(r.dpointer->privatevalue), 10) + 3);
   os << mpq_get_str(buffer, 10, r.dpointer->privatevalue);
   free(buffer);
   return os;
}



/// equality operator
bool operator==(const Rational& r, const Rational& s)
{
   return (mpq_equal(r.dpointer->privatevalue, s.dpointer->privatevalue) != 0);
}



/// inequality operator
bool operator!=(const Rational& r, const Rational& s)
{
   return (mpq_equal(r.dpointer->privatevalue, s.dpointer->privatevalue) == 0);
}



/// less than operator
bool operator<(const Rational& r, const Rational& s)
{
   return (mpq_cmp(r.dpointer->privatevalue, s.dpointer->privatevalue) < 0);
}



/// less than or equal to operator
bool operator<=(const Rational& r, const Rational& s)
{
   return (mpq_cmp(r.dpointer->privatevalue, s.dpointer->privatevalue) <= 0);
}



/// greater than operator
bool operator>(const Rational& r, const Rational& s)
{
   return (mpq_cmp(r.dpointer->privatevalue, s.dpointer->privatevalue) > 0);
}



/// greater than or equal to operator
bool operator>=(const Rational& r, const Rational& s)
{
   return (mpq_cmp(r.dpointer->privatevalue, s.dpointer->privatevalue) >= 0);
}



/// equality operator for Rational and double
bool operator==(const Rational& r, const double& s)
{
   bool res;
   mpq_t exactDouble;
   mpq_init(exactDouble);
   mpq_set_d(exactDouble, s);
   res = (mpq_equal(r.dpointer->privatevalue, exactDouble) != 0);
   mpq_clear(exactDouble);
   return res;
}



/// inequality operator for Rational and double
bool operator!=(const Rational& r, const double& s)
{
   bool res;
   mpq_t exactDouble;
   mpq_init(exactDouble);
   mpq_set_d(exactDouble, s);
   res = (mpq_equal(r.dpointer->privatevalue, exactDouble) == 0);
   mpq_clear(exactDouble);
   return res;
}



/// less than operator for Rational and double
bool operator<(const Rational& r, const double& s)
{
   bool res;
   mpq_t exactDouble;
   mpq_init(exactDouble);
   mpq_set_d(exactDouble, s);
   res = (mpq_cmp(r.dpointer->privatevalue, exactDouble) < 0);
   mpq_clear(exactDouble);
   return res;
}



/// less than or equal to operator for Rational and double
bool operator<=(const Rational& r, const double& s)
{
   bool res;
   mpq_t exactDouble;
   mpq_init(exactDouble);
   mpq_set_d(exactDouble, s);
   res = (mpq_cmp(r.dpointer->privatevalue, exactDouble) <= 0);
   mpq_clear(exactDouble);
   return res;
}



/// greater than operator for Rational and double
bool operator>(const Rational& r, const double& s)
{
   bool res;
   mpq_t exactDouble;
   mpq_init(exactDouble);
   mpq_set_d(exactDouble, s);
   res = (mpq_cmp(r.dpointer->privatevalue, exactDouble) > 0);
   mpq_clear(exactDouble);
   return res;
}



/// greater than or equal to operator for Rational and double
bool operator>=(const Rational& r, const double& s)
{
   bool res;
   mpq_t exactDouble;
   mpq_init(exactDouble);
   mpq_set_d(exactDouble, s);
   res = (mpq_cmp(r.dpointer->privatevalue, exactDouble) >= 0);
   mpq_clear(exactDouble);
   return res;
}



/// equality operator for double and Rational
bool operator==(const double& r, const Rational& s)
{
   bool res;
   mpq_t exactDouble;
   mpq_init(exactDouble);
   mpq_set_d(exactDouble, r);
   res = (mpq_equal(exactDouble, s.dpointer->privatevalue) != 0);
   mpq_clear(exactDouble);
   return res;
}



/// inequality operator double and Rational
bool operator!=(const double& r, const Rational& s)
{
   bool res;
   mpq_t exactDouble;
   mpq_init(exactDouble);
   mpq_set_d(exactDouble, r);
   res = (mpq_equal(exactDouble, s.dpointer->privatevalue) == 0);
   mpq_clear(exactDouble);
   return res;
}



/// less than operator double and Rational
bool operator<(const double& r, const Rational& s)
{
   bool res;
   mpq_t exactDouble;
   mpq_init(exactDouble);
   mpq_set_d(exactDouble, r);
   res = (mpq_cmp(exactDouble, s.dpointer->privatevalue) < 0);
   mpq_clear(exactDouble);
   return res;
}



/// less than or equal to operator double and Rational
bool operator<=(const double& r, const Rational& s)
{
   bool res;
   mpq_t exactDouble;
   mpq_init(exactDouble);
   mpq_set_d(exactDouble, r);
   res = (mpq_cmp(exactDouble, s.dpointer->privatevalue) <= 0);
   mpq_clear(exactDouble);
   return res;
}



/// greater than operator double and Rational
bool operator>(const double& r, const Rational& s)
{
   bool res;
   mpq_t exactDouble;
   mpq_init(exactDouble);
   mpq_set_d(exactDouble, r);
   res = (mpq_cmp(exactDouble, s.dpointer->privatevalue) > 0);
   mpq_clear(exactDouble);
   return res;
}



/// greater than or equal to operator double and Rational
bool operator>=(const double& r, const Rational& s)
{
   bool res;
   mpq_t exactDouble;
   mpq_init(exactDouble);
   mpq_set_d(exactDouble, r);
   res = (mpq_cmp(exactDouble, s.dpointer->privatevalue) >= 0);
   mpq_clear(exactDouble);
   return res;
}



/// addition operator for double and Rational
Rational operator+(const double& d, const Rational& r)
{
   Rational retval(d);
   mpq_add(retval.dpointer->privatevalue, retval.dpointer->privatevalue, r.dpointer->privatevalue);
   return retval;
}



/// subtraction operator for double and Rational
Rational operator-(const double& d, const Rational& r)
{
   Rational retval(d);
   mpq_sub(retval.dpointer->privatevalue, retval.dpointer->privatevalue, r.dpointer->privatevalue);
   return retval;
}



/// multiplication operator for double and Rational
Rational operator*(const double& d, const Rational& r)
{
   Rational retval(d);
   mpq_mul(retval.dpointer->privatevalue, retval.dpointer->privatevalue, r.dpointer->privatevalue);
   return retval;
}



/// division operator for double and Rational
Rational operator/(const double& d, const Rational& r)
{
   Rational retval(d);
   mpq_div(retval.dpointer->privatevalue, retval.dpointer->privatevalue, r.dpointer->privatevalue);
   return retval;
}



/// Absolute.
Rational abs(const Rational& r)
{
   Rational res;
   mpq_abs(res.dpointer->privatevalue, r.dpointer->privatevalue);
   return res;
}



/// Sign function; returns 1 if r > 0, 0 if r = 0, and -1 if r < 0.
int sign(const Rational& r)
{
      return mpq_sgn(r.dpointer->privatevalue);
}



/// Negation.
Rational operator-(const Rational& r)
{
   Rational res;
   mpq_neg(res.dpointer->privatevalue, r.dpointer->privatevalue);
   return res;
}



#elif defined(SOPLEX_WITH_GMPXX)

/// Defines the "Pimpl"-class Private
class Rational::Private : public mpq_class
{
public:

   /// default constructor
   Private()
      : mpq_class()
   {
   }

   /// copy constructor
   Private(const Private& p)
      : mpq_class()
   {
      *this = p;
   }

   /// constructor from long double
   Private(const long double& r)
      : mpq_class(double(r))
   {
   }

   /// constructor from double
   Private(const double& r)
      : mpq_class(r)
   {
   }

   /// constructor from int
   Private(const int& i)
      : mpq_class(i)
   {
   }

   /// constructor from mpq_class
   Private(const mpq_class& q)
      : mpq_class(q)
   {
   }

   /// constructor from mpq_t
   Private(const mpq_t& q)
      : mpq_class(q)
   {
   }

   /* destructor -- could be directly called but it is an internal GMP function and subject to change in further releases.
   ~Private()
   {
      ~__gmp_expr();
   }
   */
};



/// default constructor
Rational::Rational()
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private();
}



/// copy constructor
Rational::Rational(const Rational& r)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(*(r.dpointer));
}



/// copy assignment operator
Rational& Rational::operator=(const Rational &r)
{
   *(this->dpointer) = *(r.dpointer);
   return *this;
}



/// constructor from long double
Rational::Rational(const long double& r)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(r);
}



/// constructor from double
Rational::Rational(const double& r)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(r);
}



/// constructor from int
Rational::Rational(const int& i)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(i);
}


/// constructor from mpq_class
Rational::Rational(const mpq_class& q)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(q);
}



/// constructor from mpq_t
Rational::Rational(const mpq_t& q)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(q);
}



/// destructor
Rational::~Rational()
{
   dpointer->~Private();
   spx_free(dpointer);
}



/// typecasts Rational to double (allows only explicit typecast)
Rational::operator double() const
{
   return this->dpointer->get_d();
}



/// typecasts Rational to long double (allows only explicit typecast)
Rational::operator long double() const
{
   return (long double)this->dpointer->get_d();
}



/// addition operator
Rational Rational::operator+(const Rational& r) const
{
   Rational retval = *this;
   *(retval.dpointer) += *(r.dpointer);
   return retval;
}



/// addition assignment operator
Rational Rational::operator+=(const Rational& r) const
{
   *(this->dpointer) += *(r.dpointer);
   return *this;
}



/// addition operator for doubles
Rational Rational::operator+(const double& d) const
{
   Rational retval = *this;
   *(retval.dpointer) += d;
   return retval;
}



/// addition assignment operator for doubles
Rational Rational::operator+=(const double& d) const
{
   *(this->dpointer) += d;
   return *this;
}



/// subtraction operator
Rational Rational::operator-(const Rational& r) const
{
   Rational retval = *this;
   *(retval.dpointer) -= *(r.dpointer);
   return retval;
}



/// subtraction assignment operator
Rational Rational::operator-=(const Rational& r) const
{
   *(this->dpointer) -= *(r.dpointer);
   return *this;
}



/// subtraction operator for doubles
Rational Rational::operator-(const double& d) const
{
   Rational retval = *this;
   *(retval.dpointer) -= d;
   return retval;
}



/// subtraction assignment operator for doubles
Rational Rational::operator-=(const double& d) const
{
   *(this->dpointer) -= d;
   return *this;
}



/// multiplication operator
Rational Rational::operator*(const Rational& r) const
{
   Rational retval = *this;
   *(retval.dpointer) *= *(r.dpointer);
   return retval;
}



/// multiplication assignment operator
Rational Rational::operator*=(const Rational& r) const
{
   *(this->dpointer) *= *(r.dpointer);
   return *this;
}



/// multiplication operator for doubles
Rational Rational::operator*(const double& d) const
{
   Rational retval = *this;
   *(retval.dpointer) *= d;
   return retval;
}



/// multiplication assignment operator for doubles
Rational Rational::operator*=(const double& d) const
{
   *(this->dpointer) *= d;
   return *this;
}



/// division operator
Rational Rational::operator/(const Rational& r) const
{
   Rational retval = *this;
   *(retval.dpointer) /= *(r.dpointer);
   return retval;
}



/// division assignment operator
Rational Rational::operator/=(const Rational& r) const
{
   *(this->dpointer) /= *(r.dpointer);
   return *this;
}



/// division operator for doubles
Rational Rational::operator/(const double& d) const
{
   Rational retval = *this;
   *(retval.dpointer) /= d;
   return retval;
}



/// division assignment operator for doubles
Rational Rational::operator/=(const double& d) const
{
   *(this->dpointer) /= d;
   return *this;
}



/// checks if d is the closest possible double
bool Rational::isNextTo(const double& d)
{
   // get intervall [a,b] of doubles that the Rational is in
   double x = double(*this);
   double a;
   double b;

   if( x < *this )
   {
      a = x;
      b = nextafter(a, infinity);
   }
   else
   {
      b = x;
      a = nextafter(b, -infinity);
   }

   // check if d equals the closer end of the intervall
   bool result = (abs(*this - a) < abs(*this - b))
      ? (d == a)
      : (d == b);

   return result;
}



/// checks if d is exactly equal to the Rational and if not, if it is one of the two adjacent doubles
bool Rational::isAdjacentTo(const double& d)
{
   double x = double(*this);
   double a;
   double b;
   mpq_t tmp;

   mpq_init(tmp);
   mpq_set_d(tmp, x);

   int cmp = mpq_cmp(tmp, this->dpointer->get_mpq_t());
   mpq_clear(tmp);

   // the rounded value is smaller than the rational value
   if( cmp < 0 )
   {
      a = x;
      b = nextafter(a, infinity);
   }
   // the rounded value is larger than the rational value
   else if( cmp > 0 )
   {
      b = x;
      a = nextafter(b, -infinity);
   }
   // the rational is representable in double precision
   else
      return (x == d);

   return ((a == d) || (b == d));
}



/// returns precision of Rational implementation, i.e., number of bits used to store Rational numbers (INT_MAX if exact)
int Rational::precision()
{
   return INT_MAX;
}



#define MAX_STR_LEN 10000
/// read Rational from string
bool Rational::readString(const char* s)
{
   assert(s != 0);
   assert(strlen(s) <= MAX_STR_LEN);
   Rational value;
   const char* pos;

   // if there is a slash or there is no dot and exponent (i.e. we
   // have an integer), we may simply call GMP's string reader
   if( strchr(s, '/') != 0 || strpbrk(s, ".eE") == 0 )
   {
      pos = (*s == '+') ? s + 1 : s;
      if( value.dpointer->set_str(pos, 10) == 0 )
      {
         value.dpointer->canonicalize();
         *(this->dpointer) = *(value.dpointer);
         return true;
      }
      else
         return false;
   }

   // otherwise we analyze the string
#ifndef NDEBUG
   bool has_exponent = false;
   bool has_dot = false;
#endif
   bool has_digits = false;
   bool has_emptyexponent = false;
   long int exponent = 0;
   long int decshift = 0;
   mpz_class shiftpower;
   char* t;
   char tmp[MAX_STR_LEN];

   pos = s;

   // 1. sign
   if( (*pos == '+') || (*pos == '-') )
      pos++;

   // 2. Digits before the decimal dot
   while( (*pos >= '0') && (*pos <= '9') )
   {
      has_digits = true;
      pos++;
   }

   // 3. Decimal dot
   if( *pos == '.' )
   {
#ifndef NDEBUG
      has_dot = true;
#endif
      pos++;

      // 4. If there was a dot, possible digit behind it
      while( (*pos >= '0') && (*pos <= '9') )
      {
         has_digits = true;
         pos++;
      }
   }

   // 5. Exponent
   if( tolower(*pos) == 'e' )
   {
#ifndef NDEBUG
      has_exponent = true;
#endif
      has_emptyexponent = true;
      pos++;

      // 6. Exponent sign
      if( (*pos == '+') || (*pos == '-') )
         pos++;

      // 7. Exponent digits
      while( (*pos >= '0') && (*pos <= '9') )
      {
         has_emptyexponent = false;
         pos++;
      }
   }

   if( has_emptyexponent || !has_digits )
      return false;

   assert(has_exponent || has_dot);

   //read up to dot recording digits
   t = tmp;
   pos = s;

   if( *pos == '+' )
      pos++;

   while( ((*pos >= '0') && (*pos <= '9') ) || *pos == '+' || *pos == '-'  )
   {
      *t++ = *pos;
      pos++;
   }
   //record digits after dot, recording positions
   decshift = 0;
   if( *pos == '.' )
   {
      assert(has_dot);
      pos++;
      while( (*pos >= '0') && (*pos <= '9') )
      {
         *t++ = *pos;
         decshift++;
         pos++;
      }
   }
   *t = '\0';

   if( value.dpointer->set_str(tmp, 10) != 0)
      return false;
   value.dpointer->canonicalize();

   //record exponent and update final result
   exponent = -decshift;
   if( tolower(*pos) == 'e' )
   {
      pos++;
      assert(has_exponent);
      for( t = tmp; *pos != '\0'; pos++ )
         *t++ = *pos;
      *t = '\0';
      exponent += atol(tmp);
   }
   if( exponent > 0 )
   {
      mpz_ui_pow_ui(shiftpower.get_mpz_t(), 10, exponent);
      *(value.dpointer) *= shiftpower;
   }
   else if( exponent < 0 )
   {
      mpz_ui_pow_ui(shiftpower.get_mpz_t(), 10, -exponent);
      *(value.dpointer) /= shiftpower;
   }

   value.dpointer->canonicalize();
   *(this->dpointer) = *(value.dpointer);
   return true;
}



/// convert rational number to string
std::string rationalToString(const Rational& r, const bool asfloat)
{
   std::stringstream sstream;

   sstream << r;
   if( !asfloat )
      return sstream.str();
   else
   {
      sstream.str("");
      sstream << std::setprecision(20) << double(r);
      return sstream.str();
   }

   ///@todo the following code creates a string of the form 1.333...e-06; this would be nice for human readable output,
   ///      because it indicates that precision may have been lost, however, this cannot be parsed by our awk evaluation
   ///      scripts; we should think about how this could be useful
#if 0
   sstream.str("");
   sstream << mpf_class(r);
   std::string origstring = sstream.str();
   sstream.str("");

   size_t epos = origstring.find("e");

   sstream << origstring.substr(0, epos);

   if( origstring.find(".") == std::string::npos )
      sstream << ".0";

   sstream << "...";

   if( epos != std::string::npos )
      sstream << origstring.substr(origstring.find("e"));

   return sstream.str();
#endif
}



/// read Rational from string
bool readStringRational(const char* s, Rational& value)
{
   assert(s != 0);
   assert(strlen(s) <= MAX_STR_LEN);
   const char* pos;

   // if there is a slash or there is no dot and exponent (i.e. we
   // have an integer), we may simply call GMP's string reader
   if( strchr(s, '/') != 0 || strpbrk(s, ".eE") == 0 )
   {
      pos = (*s == '+') ? s + 1 : s;
      if( value.dpointer->set_str(pos, 10) == 0 )
      {
         value.dpointer->canonicalize();
         return true;
      }
      else
         return false;
   }

   // otherwise we analyze the string
#ifndef NDEBUG
   bool has_exponent = false;
   bool has_dot = false;
#endif
   bool has_digits = false;
   bool has_emptyexponent = false;
   long int exponent = 0;
   long int decshift = 0;
   mpz_class shiftpower;
   char* t;
   char tmp[MAX_STR_LEN];

   pos = s;

   // 1. sign
   if( (*pos == '+') || (*pos == '-') )
      pos++;

   // 2. Digits before the decimal dot
   while( (*pos >= '0') && (*pos <= '9') )
   {
      has_digits = true;
      pos++;
   }

   // 3. Decimal dot
   if( *pos == '.' )
   {
#ifndef NDEBUG
      has_dot = true;
#endif
      pos++;

      // 4. If there was a dot, possible digit behind it
      while( (*pos >= '0') && (*pos <= '9') )
      {
         has_digits = true;
         pos++;
      }
   }

   // 5. Exponent
   if( tolower(*pos) == 'e' )
   {
#ifndef NDEBUG
      has_exponent = true;
#endif
      has_emptyexponent = true;
      pos++;

      // 6. Exponent sign
      if( (*pos == '+') || (*pos == '-') )
         pos++;

      // 7. Exponent digits
      while( (*pos >= '0') && (*pos <= '9') )
      {
         has_emptyexponent = false;
         pos++;
      }
   }

   if( has_emptyexponent || !has_digits )
      return false;

   assert(has_exponent || has_dot);

   // read up to dot recording digits
   t = tmp;
   pos = s;

   if( *pos == '+' )
      pos++;

   while( ((*pos >= '0') && (*pos <= '9') ) || *pos == '+' || *pos == '-'  )
   {
      *t++ = *pos;
      pos++;
   }
   // record digits after dot, recording positions
   decshift = 0;
   if( *pos == '.' )
   {
      assert(has_dot);
      pos++;
      while( (*pos >= '0') && (*pos <= '9') )
      {
         *t++ = *pos;
         decshift++;
         pos++;
      }
   }
   *t = '\0';

   if( value.dpointer->set_str(tmp, 10) != 0)
      return false;

   value.dpointer->canonicalize();

   // record exponent and update final result
   exponent = -decshift;
   if( tolower(*pos) == 'e' )
   {
      pos++;
      assert(has_exponent);
      for( t = tmp; *pos != '\0'; pos++ )
         *t++ = *pos;
      *t = '\0';
      exponent += atol(tmp);
   }
   if( exponent > 0 )
   {
      mpz_ui_pow_ui(shiftpower.get_mpz_t(), 10, exponent);
      *(value.dpointer) *= shiftpower;
   }
   else if( exponent < 0 )
   {
      mpz_ui_pow_ui(shiftpower.get_mpz_t(), 10, -exponent);
      *(value.dpointer) /= shiftpower;
   }

   value.dpointer->canonicalize();
   return true;
}



/// print Rational
std::ostream& operator<<(std::ostream& os, const Rational& r)
{
   char* buffer;
   buffer = (char*) malloc (mpz_sizeinbase(mpq_numref(r.dpointer->get_mpq_t()), 10) + mpz_sizeinbase(mpq_denref(r.dpointer->get_mpq_t()), 10) + 3);
   os << mpq_get_str(buffer, 10, r.dpointer->get_mpq_t());
   free(buffer);
   return os;
}



/// equality operator
bool operator==(const Rational& r, const Rational& s)
{
   return (*(r.dpointer) == *(s.dpointer));
}



/// inequality operator
bool operator!=(const Rational& r, const Rational& s)
{
   return (*(r.dpointer) != *(s.dpointer));
}



/// less than operator
bool operator<(const Rational& r, const Rational& s)
{
   return (*(r.dpointer) < *(s.dpointer));
}



/// less than or equal to operator
bool operator<=(const Rational& r, const Rational& s)
{
   return (*(r.dpointer) <= *(s.dpointer));
}



/// greater than operator
bool operator>(const Rational& r, const Rational& s)
{
   return (*(r.dpointer) > *(s.dpointer));
}



/// greater than or equal to operator
bool operator>=(const Rational& r, const Rational& s)
{
   return (*(r.dpointer) >= *(s.dpointer));
}



/// equality operator for Rational and double
bool operator==(const Rational& r, const double& s)
{
   return (*(r.dpointer) == s);
}


/// inequality operator for Rational and double
bool operator!=(const Rational& r, const double& s)
{
   return (*(r.dpointer) != s);
}



/// less than operator for Rational and double
bool operator<(const Rational& r, const double& s)
{
   return (*(r.dpointer) < s);
}



/// less than or equal to operator for Rational and double
bool operator<=(const Rational& r, const double& s)
{
   return (*(r.dpointer) <= s);
}



/// greater than operator for Rational and double
bool operator>(const Rational& r, const double& s)
{
   return (*(r.dpointer) > s);
}



/// greater than or equal to operator for Rational and double
bool operator>=(const Rational& r, const double& s)
{
   return (*(r.dpointer) >= s);
}



/// equality operator for double and Rational
bool operator==(const double& r, const Rational& s)
{
   return (r == *(s.dpointer));
}



/// inequality operator double and Rational
bool operator!=(const double& r, const Rational& s)
{
   return (r != *(s.dpointer));
}



/// less than operator double and Rational
bool operator<(const double& r, const Rational& s)
{
   return (r < *(s.dpointer));
}



/// less than or equal to operator double and Rational
bool operator<=(const double& r, const Rational& s)
{
   return (r <= *(s.dpointer));
}



/// greater than operator double and Rational
bool operator>(const double& r, const Rational& s)
{
   return (r > *(s.dpointer));
}



/// greater than or equal to operator double and Rational
bool operator>=(const double& r, const Rational& s)
{
   return (r >= *(s.dpointer));
}



/// addition operator for double and Rational
Rational operator+(const double& d, const Rational& r)
{
   Rational retval(d);
   *(retval.dpointer) += *(r.dpointer);
   return retval;
}



/// subtraction operator for double and Rational
Rational operator-(const double& d, const Rational& r)
{
   Rational retval(d);
   *(retval.dpointer) -= *(r.dpointer);
   return retval;
}



/// multiplication operator for double and Rational
Rational operator*(const double& d, const Rational& r)
{
   Rational retval(d);
   *(retval.dpointer) *= *(r.dpointer);
   return retval;
}



/// division operator for double and Rational
Rational operator/(const double& d, const Rational& r)
{
   Rational retval(d);
   *(retval.dpointer) /= *(r.dpointer);
   return retval;
}



/// Absolute.
Rational abs(const Rational& r)
{
   Rational res = r;

   if( *(res.dpointer) < 0 )
      *(res.dpointer) *= -1;

   return res;
}



/// Sign function; returns 1 if r > 0, 0 if r = 0, and -1 if r < 0.
int sign(const Rational& r)
{
      return sgn(*(r.dpointer));
}



/// Negation
Rational operator-(const Rational& r)
{
   Rational res = r;
   *(res.dpointer) *= -1;
   return res;
}



#else



/// Defines the "Pimpl"-class Private
class Rational::Private
{

public:

   /// value
   long double privatevalue;

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

   /// copy assignment operator
   Private& operator=(const Private& p)
   {
      this->privatevalue = p.privatevalue;
      return *this;
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

};



/// default constructor
Rational::Rational()
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private();
}



/// copy constructor
Rational::Rational(const Rational& r)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(*(r.dpointer));
}



/// copy assignment operator
Rational& Rational::operator=(const Rational &r)
{
   *(this->dpointer) = *(r.dpointer);
   return *this;
}



/// constructor from long double
Rational::Rational(const long double& r)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(r);
}



/// constructor from double
Rational::Rational(const double& r)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(r);
}



/// constructor from int
Rational::Rational(const int& i)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(i);
}



/// destructor
Rational::~Rational()
{
   spx_free(dpointer);
}



/// typecasts Rational to double (allows only explicit typecast)
Rational::operator double() const
{
   return (double)this->dpointer->privatevalue;
}



/// typecasts Rational to long double (allows only explicit typecast)
Rational::operator long double() const

{
   return this->dpointer->privatevalue;
}



/// addition operator
Rational Rational::operator+(const Rational& r) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue += r.dpointer->privatevalue;
   return retval;
}



/// addition assignment operator
Rational Rational::operator+=(const Rational& r) const
{
   this->dpointer->privatevalue += r.dpointer->privatevalue;
   return *this;
}



/// addition operator for doubles
Rational Rational::operator+(const double& d) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue += d;
   return retval;
}



/// addition assignment operator for doubles
Rational Rational::operator+=(const double& d) const
{
   this->dpointer->privatevalue += d;
   return *this;
}



/// subtraction operator
Rational Rational::operator-(const Rational& r) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue -= r.dpointer->privatevalue;
   return retval;
}



/// subtraction assignment operator
Rational Rational::operator-=(const Rational& r) const
{
   this->dpointer->privatevalue -= r.dpointer->privatevalue;
   return *this;
}



/// subtraction operator for doubles
Rational Rational::operator-(const double& d) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue -= d;
   return retval;
}



/// subtraction assignment operator for doubles
Rational Rational::operator-=(const double& d) const
{
   this->dpointer->privatevalue -= d;
   return *this;
}



/// multiplication operator
Rational Rational::operator*(const Rational& r) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue *= r.dpointer->privatevalue;
   return retval;
}



/// multiplication assignment operator
Rational Rational::operator*=(const Rational& r) const
{
   this->dpointer->privatevalue *= r.dpointer->privatevalue;
   return *this;
}



/// multiplication operator for doubles
Rational Rational::operator*(const double& d) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue *= d;
   return retval;
}



/// multiplication assignment operator for doubles
Rational Rational::operator*=(const double& d) const
{
   this->dpointer->privatevalue *= d;
   return *this;
}



/// division operator
Rational Rational::operator/(const Rational& r) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue /= r.dpointer->privatevalue;
   return retval;
}



/// division assignment operator
Rational Rational::operator/=(const Rational& r) const
{
   this->dpointer->privatevalue /= r.dpointer->privatevalue;
   return *this;
}



/// division operator for doubles
Rational Rational::operator/(const double& d) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue /= d;
   return retval;
}



/// division assignment operator for doubles
Rational Rational::operator/=(const double& d) const
{
   this->dpointer->privatevalue /= d;
   return *this;
}



/// checks if d is the closest possible double
bool Rational::isNextTo(const double& d)
{
   return (this->dpointer->privatevalue == d);
}



/// checks if d is exactly equal to the Rational and if not, if it is one of the two adjacent doubles
bool Rational::isAdjacentTo(const double& d)
{
   return (this->dpointer->privatevalue == d);
}



/// returns precision of Rational implementation, i.e., number of bits used to store Rational numbers (INT_MAX if exact)
int Rational::precision()
{
   return sizeof(long double);
}



/// read Rational from string
bool Rational::readString(const char* s)
{
   return (sscanf(s, "%Lf", &this->dpointer->privatevalue) == 1 );
}




/// convert rational number to string
std::string rationalToString(const Rational& r, const bool asfloat)
{
   std::stringstream sstream;
   sstream << r;
   return sstream.str();
}



/// read Rational from string
bool readStringRational(const char* s, Rational& value)
{
   return (sscanf(s, "%Lf", &value.dpointer->privatevalue) == 1);
}



/// print Rational
std::ostream& operator<<(std::ostream& os, const Rational& r)
{
   os << r.dpointer->privatevalue;
   return os;
}



/// equality operator
bool operator==(const Rational& r, const Rational& s)
{
   return (r.dpointer->privatevalue == s.dpointer->privatevalue);
}



/// inequality operator
bool operator!=(const Rational& r, const Rational& s)
{
   return (r.dpointer->privatevalue != s.dpointer->privatevalue);
}



/// less than operator
bool operator<(const Rational& r, const Rational& s)
{
   return (r.dpointer->privatevalue < s.dpointer->privatevalue);
}



/// less than or equal to operator
bool operator<=(const Rational& r, const Rational& s)
{
   return (r.dpointer->privatevalue <= s.dpointer->privatevalue);
}



/// greater than operator
bool operator>(const Rational& r, const Rational& s)
{
   return (r.dpointer->privatevalue > s.dpointer->privatevalue);
}



/// greater than or equal to operator
bool operator>=(const Rational& r, const Rational& s)
{
   return (r.dpointer->privatevalue >= s.dpointer->privatevalue);
}



/// equality operator for Rational and double
bool operator==(const Rational& r, const double& s)
{
   return (r.dpointer->privatevalue > s - DEFAULT_EPS_ZERO)
      && (r.dpointer->privatevalue < s + DEFAULT_EPS_ZERO);
}



/// inequality operator for Rational and double
bool operator!=(const Rational& r, const double& s)
{
   return (r.dpointer->privatevalue <= s - DEFAULT_EPS_ZERO)
      || (r.dpointer->privatevalue >= s + DEFAULT_EPS_ZERO);
}



/// less than operator for Rational and double
bool operator<(const Rational& r, const double& s)
{
   return (r.dpointer->privatevalue < s);
}



/// less than or equal to operator for Rational and double
bool operator<=(const Rational& r, const double& s)
{
   return (r.dpointer->privatevalue <= s);
}



/// greater than operator for Rational and double
bool operator>(const Rational& r, const double& s)
{
   return (r.dpointer->privatevalue > s);
}



/// greater than or equal to operator for Rational and double
bool operator>=(const Rational& r, const double& s)
{
   return (r.dpointer->privatevalue >= s);
}



/// equality operator for double and Rational
bool operator==(const double& r, const Rational& s)
{
   return (s.dpointer->privatevalue > r - DEFAULT_EPS_ZERO)
      && (s.dpointer->privatevalue < r + DEFAULT_EPS_ZERO);
}



/// inequality operator double and Rational
bool operator!=(const double& r, const Rational& s)
{
   return (s.dpointer->privatevalue <= r - DEFAULT_EPS_ZERO)
      || (s.dpointer->privatevalue >= r + DEFAULT_EPS_ZERO);
}



/// less than operator double and Rational
bool operator<(const double& r, const Rational& s)
{
   return (r < s.dpointer->privatevalue);
}



/// less than or equal to operator double and Rational
bool operator<=(const double& r, const Rational& s)
{
   return (r <= s.dpointer->privatevalue);
}



/// greater than operator double and Rational
bool operator>(const double& r, const Rational& s)
{
   return (r > s.dpointer->privatevalue);
}



/// greater than or equal to operator double and Rational
bool operator>=(const double& r, const Rational& s)
{
   return (r >= s.dpointer->privatevalue);
}



/// addition operator for double and Rational
Rational operator+(const double& d, const Rational& r)
{
   Rational retval(d);
   retval.dpointer->privatevalue += r.dpointer->privatevalue;
   return retval;
}



/// subtraction operator for double and Rational
Rational operator-(const double& d, const Rational& r)
{
   Rational retval(d);
   retval.dpointer->privatevalue -= r.dpointer->privatevalue;
   return retval;
}



/// multiplication operator for double and Rational
Rational operator*(const double& d, const Rational& r)
{
   Rational retval(d);
   retval.dpointer->privatevalue *= r.dpointer->privatevalue;
   return retval;
}



/// division operator for double and Rational
Rational operator/(const double& d, const Rational& r)
{
   Rational retval(d);
   retval.dpointer->privatevalue /= r.dpointer->privatevalue;
   return retval;
}



/// Absolute.
Rational abs(const Rational& r)
{
   Rational res = r;

   if( res.dpointer->privatevalue < 0 )
      res.dpointer->privatevalue *= -1;

   return res;
}



/// Sign function; returns 1 if r > 0, 0 if r = 0, and -1 if r < 0.
int sign(const Rational& r)
{
      return (r.dpointer->privatevalue > 0) - (r.dpointer->privatevalue < 0);
}



/// Negation.
Rational operator-(const Rational& r)
{
   Rational res = r;
   res.dpointer->privatevalue *= -1;
   return res;
}
#endif // SOPLEX_WITH_GMP
} // namespace soplex
#endif
