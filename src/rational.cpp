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

/**@file  rational.cpp
 * @brief Wrapper for GMP types.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <string>
#include <sstream>


#include "rational.h"
#include "spxalloc.h"

#ifdef SOPLEX_WITH_GMP
#include "gmp.h"
#include "gmpxx.h"
#endif

namespace soplex
{
#ifdef SOPLEX_WITH_GMP
#define RationalIsExact() (true)

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

#if 0

   Private(const mpq_class& q)
      : mpq_class(q)
   {
   }

   Private(const mpq_t& q)
      : mpq_class(q)
   {
   }

#endif

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

///copy assignment constructor
Rational& Rational::operator=(const Rational &r)
{
   *dpointer = *(r.dpointer);
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

#if 0

Rational::Rational(const mpq_class& q)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(q);
}

Rational::Rational(const mpq_t& q)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(q);
}

#endif

/// typecasts Rational to double (allows only explicit typecast)
Rational::operator double() const
{
   return this->dpointer->get_d();
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
}

/// subtraction operator
Rational Rational::operator-(const Rational& r) const
{
   Rational retval = *this;
   *(retval.dpointer) -= *(r.dpointer);
   return retval;
}

// subtraction assignment operator
Rational Rational::operator-=(const Rational& r) const
{
   *(this->dpointer) -= *(r.dpointer);
}

/// subtraction operator for doubles
Rational Rational::operator-(const double& d) const
{
   Rational retval = *this;
   *(retval.dpointer) -= d;
   return retval;
}

// subtraction assignment operator for doubles
Rational Rational::operator-=(const double& d) const
{
   *(this->dpointer) -= d;
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
   bool has_digits = false;
   bool has_exponent = false;
   bool has_dot = false;
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
      has_dot = true;
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
      has_exponent = true;
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

   assert( has_exponent || has_dot);

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
      sstream << std::setprecision(20) << mpf_class(*(r.dpointer), 128);
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
   bool has_digits = false;
   bool has_exponent = false;
   bool has_dot = false;
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
      has_dot = true;
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
      has_exponent = true;
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

   assert( has_exponent || has_dot);

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

/// Absolute.
Rational abs(const Rational& r)
{
   Rational res = r;

   if( *(res.dpointer) < 0 )
      *(res.dpointer) *= -1;

   return res;
}

/// print Rational
std::ostream& operator<<(std::ostream& os, const Rational& q)
{
   os << mpq_class(*(q.dpointer));
   return os;
}

/// Negation.
Rational operator-(const Rational& q)
{
   Rational res = q;
   *(res.dpointer) *= -1;
   return res;
}

/// equality operator
bool operator==(const Rational& r, const Rational& s)
{
   bool res = (*(r.dpointer) == *(s.dpointer));
   return res;
}

/// inequality operator
bool operator!=(const Rational& r, const Rational& s)
{
   bool res = (*(r.dpointer) != *(s.dpointer));
   return res;
}

/// less than operator
bool operator<(const Rational& r, const Rational& s)
{
   bool res = (*(r.dpointer) < *(s.dpointer));
   return res;
}

/// less than or equal to operator
bool operator<=(const Rational& r, const Rational& s)
{
   bool res = (*(r.dpointer) <= *(s.dpointer));
   return res;
}

/// greater than operator
bool operator>(const Rational& r, const Rational& s)
{
   bool res = (*(r.dpointer) > *(s.dpointer));
   return res;
}

/// greater than or equal to operator
bool operator>=(const Rational& r, const Rational& s)
{
   bool res = (*(r.dpointer) >= *(s.dpointer));
   return res;
}



/// equality operator for Rational and double
bool operator==(const Rational& r, const double& s)
{
   bool res = (*(r.dpointer) == s);
   return res;
}

/// inequality operator for Rational and double
bool operator!=(const Rational& r, const double& s)
{
   bool res = (*(r.dpointer) != s);
   return res;
}

/// less than operator for Rational and double
bool operator<(const Rational& r, const double& s)
{
   bool res = (*(r.dpointer) < s);
   return res;
}

/// less than or equal to operator for Rational and double
bool operator<=(const Rational& r, const double& s)
{
   bool res = (*(r.dpointer) <= s);
   return res;
}

/// greater than operator for Rational and double
bool operator>(const Rational& r, const double& s)
{
   bool res = (*(r.dpointer) > s);
   return res;
}

/// greater than or equal to operator for Rational and double
bool operator>=(const Rational& r, const double& s)
{
   bool res = (*(r.dpointer) >= s);
   return res;
}



/// equality operator for double and Rational
bool operator==(const double& r, const Rational& s)
{
   bool res = (r == *(s.dpointer));
   return res;
}

/// inequality operator double and Rational
bool operator!=(const double& r, const Rational& s)
{
   bool res = (r != *(s.dpointer));
   return res;
}

/// less than operator double and Rational
bool operator<(const double& r, const Rational& s)
{
   bool res = (r < *(s.dpointer));
   return res;
}

/// less than or equal to operator double and Rational
bool operator<=(const double& r, const Rational& s)
{
   bool res = (r <= *(s.dpointer));
   return res;
}

/// greater than operator double and Rational
bool operator>(const double& r, const Rational& s)
{
   bool res = (r > *(s.dpointer));
   return res;
}

/// greater than or equal to operator double and Rational
bool operator>=(const double& r, const Rational& s)
{
   bool res = (r >= *(s.dpointer));
   return res;
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


#else
#define RationalIsExact() (false)

/// Defines the "Pimpl"-class Private
class Rational::Private
{

public:

   /// value
   double privatevalue;

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

#if 0

   Private(const mpq_class& q)
      : mpq_class(q)
   {
   }

   Private(const mpq_t& q)
      : mpq_class(q)
   {
   }

#endif

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

///copy assignment constructor
Rational& Rational::operator=(const Rational &r)
{
   *dpointer = *(r.dpointer);
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

#if 0

Rational::Rational(const mpq_class& q)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(q);
}

Rational::Rational(const mpq_t& q)
{
   dpointer = 0;
   spx_alloc(dpointer);
   dpointer = new (dpointer) Private(q);
}

#endif

/// typecasts Rational to double (allows only explicit typecast)

Rational::operator double() const

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
}

/// subtraction operator
Rational Rational::operator-(const Rational& r) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue -= r.dpointer->privatevalue;
   return retval;
}

// subtraction assignment operator
Rational Rational::operator-=(const Rational& r) const
{
   this->dpointer->privatevalue -= r.dpointer->privatevalue;
}

/// subtraction operator for doubles
Rational Rational::operator-(const double& d) const
{
   Rational retval = *this;
   retval.dpointer->privatevalue -= d;
   return retval;
}

// subtraction assignment operator for doubles
Rational Rational::operator-=(const double& d) const
{
   this->dpointer->privatevalue -= d;
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
}

/// read Rational from string
bool Rational::readString(const char* s)
{
   return (sscanf(s, "%" REAL_FORMAT, &this->dpointer->privatevalue) == 1 );
}


/// convert rational number to string
std::string rationalToString(const Rational& r, const bool asfloat)
{
   std::stringstream sstream;
   sstream << r;
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


/// read Rational from string
bool readStringRational(const char* s, Rational& value)
{
   return (sscanf(s, "%" REAL_FORMAT, &value.dpointer->privatevalue) == 1);
}

/// Absolute.
Rational abs(const Rational& r)
{
   Rational res = r;

   if( res.dpointer->privatevalue < 0 )
      res.dpointer->privatevalue *= -1;

   return res;
}

/// print Rational
std::ostream& operator<<(std::ostream& os, const Rational& q)
{
   os << q.dpointer->privatevalue;
   return os;
}

/// Negation.
Rational operator-(const Rational& q)
{
   Rational res = q;
   res.dpointer->privatevalue *= -1;
   return res;
}

/// equality operator
bool operator==(const Rational& r, const Rational& s)
{
   bool res = (r.dpointer->privatevalue == s.dpointer->privatevalue);
   return res;
}

/// inequality operator
bool operator!=(const Rational& r, const Rational& s)
{
   bool res = (r.dpointer->privatevalue != s.dpointer->privatevalue);
   return res;
}

/// less than operator
bool operator<(const Rational& r, const Rational& s)
{
   bool res = (r.dpointer->privatevalue < s.dpointer->privatevalue);
   return res;
}

/// less than or equal to operator
bool operator<=(const Rational& r, const Rational& s)
{
   bool res = (r.dpointer->privatevalue <= s.dpointer->privatevalue);
   return res;
}

/// greater than operator
bool operator>(const Rational& r, const Rational& s)
{
   bool res = (r.dpointer->privatevalue > s.dpointer->privatevalue);
   return res;
}

/// greater than or equal to operator
bool operator>=(const Rational& r, const Rational& s)
{
   bool res = (r.dpointer->privatevalue >= s.dpointer->privatevalue);
   return res;
}



/// equality operator for Rational and double
bool operator==(const Rational& r, const double& s)
{
   bool res = (r.dpointer->privatevalue == s);
   return res;
}

/// inequality operator for Rational and double
bool operator!=(const Rational& r, const double& s)
{
   bool res = (r.dpointer->privatevalue != s);
   return res;
}

/// less than operator for Rational and double
bool operator<(const Rational& r, const double& s)
{
   bool res = (r.dpointer->privatevalue < s);
   return res;
}

/// less than or equal to operator for Rational and double
bool operator<=(const Rational& r, const double& s)
{
   bool res = (r.dpointer->privatevalue <= s);
   return res;
}

/// greater than operator for Rational and double
bool operator>(const Rational& r, const double& s)
{
   bool res = (r.dpointer->privatevalue > s);
   return res;
}

/// greater than or equal to operator for Rational and double
bool operator>=(const Rational& r, const double& s)
{
   bool res = (r.dpointer->privatevalue >= s);
   return res;
}



/// equality operator for double and Rational
bool operator==(const double& r, const Rational& s)
{
   bool res = (r == s.dpointer->privatevalue);
   return res;
}

/// inequality operator double and Rational
bool operator!=(const double& r, const Rational& s)
{
   bool res = (r != s.dpointer->privatevalue);
   return res;
}

/// less than operator double and Rational
bool operator<(const double& r, const Rational& s)
{
   bool res = (r < s.dpointer->privatevalue);
   return res;
}

/// less than or equal to operator double and Rational
bool operator<=(const double& r, const Rational& s)
{
   bool res = (r <= s.dpointer->privatevalue);
   return res;
}

/// greater than operator double and Rational
bool operator>(const double& r, const Rational& s)
{
   bool res = (r > s.dpointer->privatevalue);
   return res;
}

/// greater than or equal to operator double and Rational
bool operator>=(const double& r, const Rational& s)
{
   bool res = (r >= s.dpointer->privatevalue);
   return res;
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



#endif
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
