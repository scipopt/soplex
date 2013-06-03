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

#include "rational.h"
#include <string>
#include <sstream>

namespace soplex
{
#ifdef SOPLEX_WITH_GMP

/// Negation.
Rational operator-(const Rational& q)
{
   Rational res = q;
   res *= -1;
   return res;
}

/// Division.
Rational operator/(const Rational& p, const Rational& q)
{
   Rational res = p;
   res /= q;
   return res;
}

/// print Rational
std::ostream& operator<<(std::ostream& os, const Rational& q)
{
   os << mpq_class(q);
   return os;
}

/// return as string
std::string Rational::toString(const bool asfloat)
{
   std::stringstream sstream;

   if( asfloat )
      sstream << mpf_class(*this);
   else
      sstream << mpq_class(*this);

   return sstream.str();
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
      if( value.set_str(pos, 10) == 0 )
      {
         value.canonicalize();
         *this = value;
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

   if( value.set_str(tmp, 10) != 0)
      return false;
   value.canonicalize();

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
      value *= shiftpower;
   }
   else if( exponent < 0 )
   {
      mpz_ui_pow_ui(shiftpower.get_mpz_t(), 10, -exponent);
      value /= shiftpower;
   }

   value.canonicalize();
   *this = value;
   return true;
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
