/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  spxdefines.h
 * @brief Debugging, floating point type and parameter definitions.
 *
 * In optimized code with \c NDEBUG defined, only
 * \ref soplex::SPxOut::VERB_INFO1 "VERB_INFO1",
 * \ref soplex::SPxOut::VERB_INFO2 "VERB_INFO2", and
 * \ref soplex::SPxOut::VERB_INFO3 "VERB_INFO3" are set.
 * If \c NDEBUG is not defined, the code within \#TRACE is used.
 * If \c SOPLEX_DEBUG is defined, the code within
 * \ref soplex::SPxOut::VERB_DEBUG "VERB_DEBUG" is also used.
 *
 * If \c WITH_LONG_DOUBLE is defined, all Real numbers are of type
 * long double instead of just double.
 */
#ifndef _SPXDEFINES_H_
#define _SPXDEFINES_H_

/*
 * include build configuration flags
 */
#ifndef SOPLEX_NO_CONFIG_HEADER
#include "soplex/config.h"
#endif

#ifdef SOPLEX_WITH_BOOST
#include "boost/multiprecision/number.hpp"
#ifdef SOPLEX_WITH_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

#ifdef SOPLEX_WITH_MPFR
// For multiple precision
#include <boost/multiprecision/mpfr.hpp>
#ifndef NDEBUG
#include "boost/multiprecision/debug_adaptor.hpp" // For debuging mpf numbers
#endif // NDEBUG
#endif // SOPLEX_WITH_MPFR
#ifdef SOPLEX_WITH_CPPMPF
#include <boost/serialization/nvp.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif  // SOPLEX_WITH_CPPMPF

#ifdef SOPLEX_WITH_GMP
#include <boost/multiprecision/gmp.hpp>
#else
#include <boost/multiprecision/cpp_int.hpp>
#endif

#endif

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <memory>
#include <limits>
#include <cmath>
#ifdef _MSC_VER
#include <float.h>
#endif


namespace soplex
{
// Overloaded EQ function
bool EQ(int a, int b);

#define SOPLEX_VERSION         900
#define SOPLEX_VERSION_SUB       0  ///< @deprecated Always 0
#define SOPLEX_SUBVERSION        0  ///< @deprecated Always 0
#define SOPLEX_APIVERSION       20
#define SOPLEX_COPYRIGHT       "Copyright (c) 1996-2025 Zuse Institute Berlin (ZIB)"

/*-----------------------------------------------------------------------------
 * Assertion Macros etc.
 *-----------------------------------------------------------------------------
 */

/**
   \brief Macro to turn some assertions into warnings.

   If both \c NDEBUG and \c WITH_WARNINGS are defined then the failed
   assertion is converted to a warning. In all other cases this macro is
   equivalent to assert().

   @param  prefix  Short string for grepping in source code.
   @param  expr    Expression that must be satisfied.
*/
#if defined (NDEBUG) && defined (WITH_WARNINGS)
#define SOPLEX_ASSERT_WARN( prefix, expr )                        \
   if ( !( expr ) )                                        \
      {                                                    \
         std::cerr                                         \
         << prefix                                         \
         << " failed assertion on line " << __LINE__       \
         << " in file " << __FILE__ << ": "                \
         << #expr                                          \
         << std::endl;                                     \
      }
#else // just a normal assert
#define SOPLEX_ASSERT_WARN( prefix, expr ) ( assert( expr ) )
#endif



/*-----------------------------------------------------------------------------
 * Debugging Macros etc.
 *-----------------------------------------------------------------------------
 */

/**
   Prints/Executes \p stream with verbosity level \p verbosity, resetting
   the old verbosity level afterwards.
   Usually the parameter \p stream prints something out.
   This is an internal define used by SPX_MSG_ERROR, SPX_MSG_WARNING, etc.
*/
#ifdef DISABLE_VERBOSITY
#define SOPLEX_DO_WITH_TMP_VERBOSITY( verbosity, spxout, do_something ) {}
#define SOPLEX_DO_WITH_ERR_VERBOSITY( do_something ) {}
#else
#define SOPLEX_DO_WITH_TMP_VERBOSITY( verbosity, spxout, do_something ) \
   {                                                             \
     if( &spxout != nullptr )                                    \
     {                                                           \
        if( verbosity <= spxout.getVerbosity() )                 \
        {                                                        \
           const SPxOut::Verbosity  old_verbosity = spxout.getVerbosity(); \
           spxout.setVerbosity( verbosity );                     \
           do_something;                                         \
           spxout.setVerbosity( old_verbosity );                 \
        }                                                        \
     }                                                           \
   }
#define SOPLEX_DO_WITH_ERR_VERBOSITY( do_something ) { do_something; }
#endif

/// Prints out message \p x if the verbosity level is at least SPxOut::VERB_ERROR.
#define SPX_MSG_ERROR(x)            { SOPLEX_DO_WITH_ERR_VERBOSITY( x ) }
/// Prints out message \p x if the verbosity level is at least SPxOut::VERB_WARNING.
#define SPX_MSG_WARNING(spxout, x)  { SOPLEX_DO_WITH_TMP_VERBOSITY( SPxOut::VERB_WARNING, spxout, x ) }
/// Prints out message \p x if the verbosity level is at least SPxOut::VERB_INFO1.
#define SPX_MSG_INFO1(spxout, x)    { SOPLEX_DO_WITH_TMP_VERBOSITY( SPxOut::VERB_INFO1, spxout, x ) }
/// Prints out message \p x if the verbosity level is at least SPxOut::VERB_INFO2.
#define SPX_MSG_INFO2(spxout, x)    { SOPLEX_DO_WITH_TMP_VERBOSITY( SPxOut::VERB_INFO2, spxout, x ) }
/// Prints out message \p x if the verbosity level is at least SPxOut::VERB_INFO3.
#define SPX_MSG_INFO3(spxout, x)    { SOPLEX_DO_WITH_TMP_VERBOSITY( SPxOut::VERB_INFO3, spxout, x ) }

extern bool msginconsistent(const char* name, const char* file, int line);

#define SPX_MSG_INCONSISTENT(name) msginconsistent(name, __FILE__, __LINE__)

#if defined(SOPLEX_DEBUG)
// print output in any case, regardless of _tolerances->verbose():
#define SPX_MSG_DEBUG(x) { x; }
#define SPX_DEBUG(x) { x; }
#else
#define SPX_MSG_DEBUG(x) /**/
#define SPX_DEBUG(x) /**/
#endif //!SOPLEX_DEBUG


/*-----------------------------------------------------------------------------
 * multi-thread support
 *-----------------------------------------------------------------------------
 */
// enable the user to compile without thread_local by setting USRCXXFLAGS=-DTHREADLOCAL=""
#if !defined(SOPLEX_THREADLOCAL)
#if defined(_MSC_VER) && _MSC_VER < 1900
#define SOPLEX_THREADLOCAL
#else
#define SOPLEX_THREADLOCAL thread_local
#endif
#endif

/*-----------------------------------------------------------------------------
 * Long double support, Parameters and Epsilons
 *-----------------------------------------------------------------------------
 */


#ifdef WITH_LONG_DOUBLE


typedef long double Real;

#ifndef SOPLEX_REAL
#define SOPLEX_REAL(x)  x##L
#define SOPLEX_REAL_FORMAT "Lf"
#endif
/// default allowed bound violation
#ifndef SOPLEX_DEFAULT_BND_VIOL
#define SOPLEX_DEFAULT_BND_VIOL   1e-12L
#endif
/// default allowed additive zero: 1.0 + EPS_ZERO == 1.0
#ifndef SOPLEX_DEFAULT_EPS_ZERO
#define SOPLEX_DEFAULT_EPS_ZERO   1e-28L
#endif
/// epsilon for factorization
#ifndef SOPLEX_DEFAULT_EPS_FACTOR
#define SOPLEX_DEFAULT_EPS_FACTOR 1e-30L
#endif
/// epsilon for factorization update
#ifndef SOPLEX_DEFAULT_EPS_UPDATE
#define SOPLEX_DEFAULT_EPS_UPDATE 1e-26L
#endif
#ifndef SOPLEX_DEFAULT_EPS_PIVOR
#define SOPLEX_DEFAULT_EPS_PIVOR 1e-20L
#endif
///
#define SOPLEX_DEFAULT_INFINITY   1e100L


#else

#ifdef WITH_FLOAT

typedef float Real;

#ifndef SOPLEX_REAL
#define SOPLEX_REAL(x)  x
#define SOPLEX_REAL_FORMAT "f"
#endif
/// default allowed bound violation
#ifndef SOPLEX_DEFAULT_BND_VIOL
#define SOPLEX_DEFAULT_BND_VIOL   1e-1f
#endif
/// default allowed additive zero: 1.0 + EPS_ZERO == 1.0
#ifndef SOPLEX_DEFAULT_EPS_ZERO
#define SOPLEX_DEFAULT_EPS_ZERO   1e-7f
#endif
#ifndef SOPLEX_DEFAULT_EPS_FACTOR
#define SOPLEX_DEFAULT_EPS_FACTOR 1e-7f
#endif
#ifndef SOPLEX_DEFAULT_EPS_UPDATE
#define SOPLEX_DEFAULT_EPS_UPDATE 1e-6f
#endif
#ifndef SOPLEX_DEFAULT_EPS_PIVOR
#define SOPLEX_DEFAULT_EPS_PIVOR 1e-6f
#endif
#define SOPLEX_DEFAULT_INFINITY   1e35f

#else

typedef double Real;

#ifndef SOPLEX_REAL
#define SOPLEX_REAL(x)  x
#define SOPLEX_REAL_FORMAT "lf"
#endif
/// default allowed bound violation
#ifndef SOPLEX_DEFAULT_BND_VIOL
#define SOPLEX_DEFAULT_BND_VIOL   1e-6
#endif
/// default allowed additive zero: 1.0 + EPS_ZERO == 1.0
#ifndef SOPLEX_DEFAULT_EPS_ZERO
#define SOPLEX_DEFAULT_EPS_ZERO   1e-16
#endif
#ifndef SOPLEX_DEFAULT_EPS_FACTOR
#define SOPLEX_DEFAULT_EPS_FACTOR 1e-20
#endif
#ifndef SOPLEX_DEFAULT_EPS_UPDATE
#define SOPLEX_DEFAULT_EPS_UPDATE 1e-16
#endif
#ifndef SOPLEX_DEFAULT_EPS_PIVOR
#define SOPLEX_DEFAULT_EPS_PIVOR 1e-10
#endif
#define SOPLEX_DEFAULT_INFINITY   1e100

#endif // !WITH_FLOAT
#endif // !WITH_LONG_DOUBLE

#define SOPLEX_MAX(x,y)        ((x)>(y) ? (x) : (y))
#define SOPLEX_MIN(x,y)        ((x)<(y) ? (x) : (y))

#define SPX_MAXSTRLEN       1024 /**< maximum string length in SoPlex */

SOPLEX_THREADLOCAL extern const Real infinity;

/**
 * @brief Helper template for converting a string to a numeric type.
 *
 * General template uses R(str.c_str()) which works for boost multiprecision types.
 * Specialization for Real uses std::stod/std::stold/std::stof depending on Real type.
 */
template <class R>
struct StringToNumber
{
   static R convert(const std::string& str)
   {
      return R(str.c_str());
   }
};

/**
 * @brief Helper template for computing square root of a numeric type.
 *
 * For MPFR/cpp_dec_float types, uses boost::multiprecision::sqrt.
 * For rational types (which don't support exact sqrt), converts through double.
 * For Real (double/float/long double), uses std::sqrt.
 */
template <class R>
struct SafeSqrt
{
   static R compute(const R& val)
   {
      return sqrt(val);
   }
};

/// Specialization for Real types
template <>
struct SafeSqrt<Real>
{
   static Real compute(const Real& val)
   {
      return std::sqrt(val);
   }
};

#ifdef SOPLEX_WITH_BOOST
#ifdef SOPLEX_WITH_GMP
/// Specialization for GMP rationals - convert through double since rational sqrt isn't exact
template <>
struct SafeSqrt<boost::multiprecision::number<boost::multiprecision::gmp_rational, boost::multiprecision::et_off>>
{
   static boost::multiprecision::number<boost::multiprecision::gmp_rational, boost::multiprecision::et_off>
   compute(
      const boost::multiprecision::number<boost::multiprecision::gmp_rational, boost::multiprecision::et_off>&
      val)
   {
      // Convert to double, take sqrt, convert back (approximate but avoids irrational results)
      double dval = static_cast<double>(val);
      return boost::multiprecision::number<boost::multiprecision::gmp_rational, boost::multiprecision::et_off>
             (std::sqrt(dval));
   }
};
#endif

#ifndef SOPLEX_WITH_GMP
/// Specialization for cpp_rational (when GMP is not available)
template <>
struct SafeSqrt<boost::multiprecision::cpp_rational>
{
   static boost::multiprecision::cpp_rational compute(const boost::multiprecision::cpp_rational& val)
   {
      double dval = static_cast<double>(val);
      return boost::multiprecision::cpp_rational(std::sqrt(dval));
   }
};
#endif
#endif

#ifdef WITH_LONG_DOUBLE
/// Specialization for Real = long double
template <>
struct StringToNumber<Real>
{
   static Real convert(const std::string& str)
   {
      return std::stold(str);
   }
};
#elif defined(WITH_FLOAT)
/// Specialization for Real = float
template <>
struct StringToNumber<Real>
{
   static Real convert(const std::string& str)
   {
      return std::stof(str);
   }
};
#else
/// Specialization for Real = double
template <>
struct StringToNumber<Real>
{
   static Real convert(const std::string& str)
   {
      return std::stod(str);
   }
};
#endif

/**
 * @brief Precision traits for computing default tolerance values based on numeric type.
 *
 * This template provides precision-appropriate default tolerances for arbitrary numeric types.
 * Specializations exist for Real (double/long double/float) and Boost multiprecision types.
 */
template <class R>
struct PrecisionTraits
{
   /// Default epsilon for general zero comparisons
   static R defaultEpsilon()
   {
      // Generic fallback: use Real defaults converted to R
      return R(SOPLEX_DEFAULT_EPS_ZERO);
   }

   /// Default epsilon for factorization
   static R defaultEpsilonFactorization()
   {
      return R(SOPLEX_DEFAULT_EPS_FACTOR);
   }

   /// Default epsilon for factorization update
   static R defaultEpsilonUpdate()
   {
      return R(SOPLEX_DEFAULT_EPS_UPDATE);
   }

   /// Default epsilon for pivot tolerance
   static R defaultEpsilonPivot()
   {
      return R(SOPLEX_DEFAULT_EPS_PIVOR);
   }

   /// Default feasibility/optimality tolerance
   static R defaultFeastol()
   {
      return R(SOPLEX_DEFAULT_BND_VIOL);
   }

   /// Default infinity value
   static R defaultInfinity()
   {
      return R(SOPLEX_DEFAULT_INFINITY);
   }
};

/// Specialization of PrecisionTraits for Real (double/long double/float)
template <>
struct PrecisionTraits<Real>
{
   static Real defaultEpsilon()
   {
      return SOPLEX_DEFAULT_EPS_ZERO;
   }
   static Real defaultEpsilonFactorization()
   {
      return SOPLEX_DEFAULT_EPS_FACTOR;
   }
   static Real defaultEpsilonUpdate()
   {
      return SOPLEX_DEFAULT_EPS_UPDATE;
   }
   static Real defaultEpsilonPivot()
   {
      return SOPLEX_DEFAULT_EPS_PIVOR;
   }
   static Real defaultFeastol()
   {
      return SOPLEX_DEFAULT_BND_VIOL;
   }
   static Real defaultInfinity()
   {
      return SOPLEX_DEFAULT_INFINITY;
   }
};

#ifdef SOPLEX_WITH_BOOST

/**
 * @brief PrecisionTraits specialization for Boost multiprecision MPFR types.
 *
 * Computes precision-appropriate default tolerances based on the number of
 * precision digits available. For example, a 512-bit (155 decimal digit) type
 * would use epsilon around 1e-120, while a 1024-bit type would use ~1e-280.
 */
#ifdef SOPLEX_WITH_MPFR
template <unsigned Digits, boost::multiprecision::mpfr_allocation_type AllocType, boost::multiprecision::expression_template_option ET>
struct PrecisionTraits<boost::multiprecision::number<
   boost::multiprecision::mpfr_float_backend<Digits, AllocType>, ET>>
{
   using R = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<Digits, AllocType>, ET>;

   /// Epsilon at ~80% of available precision
   static R defaultEpsilon()
   {
      // Digits is in decimal digits; use 80% of available precision for epsilon
      int expDigits = static_cast<int>(Digits * 0.8);
      R result = R(1);

      for(int i = 0; i < expDigits; ++i)
         result /= R(10);

      return result;
   }

   /// Epsilon for factorization: tighter by factor of 1e-4
   static R defaultEpsilonFactorization()
   {
      return defaultEpsilon() * R("1e-4");
   }

   /// Epsilon for factorization update: same as general epsilon
   static R defaultEpsilonUpdate()
   {
      return defaultEpsilon();
   }

   /// Epsilon for pivot: more relaxed by factor of 1e6
   static R defaultEpsilonPivot()
   {
      int expDigits = static_cast<int>(Digits * 0.6);
      R result = R(1);

      for(int i = 0; i < expDigits; ++i)
         result /= R(10);

      return result;
   }

   /// Feasibility tolerance at ~65% of available precision
   static R defaultFeastol()
   {
      int expDigits = static_cast<int>(Digits * 0.65);
      R result = R(1);

      for(int i = 0; i < expDigits; ++i)
         result /= R(10);

      return result;
   }

   /// Infinity at ~50% of available precision (to avoid overflow)
   static R defaultInfinity()
   {
      int expDigits = static_cast<int>(Digits * 0.5);
      R result = R(1);

      for(int i = 0; i < expDigits; ++i)
         result *= R(10);

      return result;
   }
};

/// Specialization for MPFR with 0 digits (variable precision at runtime)
template <boost::multiprecision::mpfr_allocation_type AllocType, boost::multiprecision::expression_template_option ET>
struct PrecisionTraits<boost::multiprecision::number<
   boost::multiprecision::mpfr_float_backend<0, AllocType>, ET>>
{
   using R = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0, AllocType>, ET>;

   static R defaultEpsilon()
   {
      // For variable precision MPFR, use the current default precision
      unsigned prec = boost::multiprecision::mpfr_float_backend<0>::default_precision();
      // Convert bits to decimal digits: digits10 ~ prec * log10(2) ~ prec * 0.301
      int decDigits = static_cast<int>(prec * 0.301 * 0.8);
      R result = R(1);

      for(int i = 0; i < decDigits; ++i)
         result /= R(10);

      return result;
   }

   static R defaultEpsilonFactorization()
   {
      return defaultEpsilon() * R("1e-4");
   }

   static R defaultEpsilonUpdate()
   {
      return defaultEpsilon();
   }

   static R defaultEpsilonPivot()
   {
      unsigned prec = boost::multiprecision::mpfr_float_backend<0>::default_precision();
      int decDigits = static_cast<int>(prec * 0.301 * 0.6);
      R result = R(1);

      for(int i = 0; i < decDigits; ++i)
         result /= R(10);

      return result;
   }

   static R defaultFeastol()
   {
      unsigned prec = boost::multiprecision::mpfr_float_backend<0>::default_precision();
      int decDigits = static_cast<int>(prec * 0.301 * 0.65);
      R result = R(1);

      for(int i = 0; i < decDigits; ++i)
         result /= R(10);

      return result;
   }

   static R defaultInfinity()
   {
      unsigned prec = boost::multiprecision::mpfr_float_backend<0>::default_precision();
      int decDigits = static_cast<int>(prec * 0.301 * 0.5);
      R result = R(1);

      for(int i = 0; i < decDigits; ++i)
         result *= R(10);

      return result;
   }
};
#endif // SOPLEX_WITH_MPFR

#ifdef SOPLEX_WITH_CPPMPF
/// Specialization for cpp_dec_float (Boost's pure C++ decimal float)
template <unsigned Digits, boost::multiprecision::expression_template_option ET>
struct PrecisionTraits<boost::multiprecision::number<
   boost::multiprecision::cpp_dec_float<Digits>, ET>>
{
   using R = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<Digits>, ET>;

   static R defaultEpsilon()
   {
      int expDigits = static_cast<int>(Digits * 0.8);
      R result = R(1);

      for(int i = 0; i < expDigits; ++i)
         result /= R(10);

      return result;
   }

   static R defaultEpsilonFactorization()
   {
      return defaultEpsilon() * R("1e-4");
   }

   static R defaultEpsilonUpdate()
   {
      return defaultEpsilon();
   }

   static R defaultEpsilonPivot()
   {
      int expDigits = static_cast<int>(Digits * 0.6);
      R result = R(1);

      for(int i = 0; i < expDigits; ++i)
         result /= R(10);

      return result;
   }

   static R defaultFeastol()
   {
      int expDigits = static_cast<int>(Digits * 0.65);
      R result = R(1);

      for(int i = 0; i < expDigits; ++i)
         result /= R(10);

      return result;
   }

   static R defaultInfinity()
   {
      int expDigits = static_cast<int>(Digits * 0.5);
      R result = R(1);

      for(int i = 0; i < expDigits; ++i)
         result *= R(10);

      return result;
   }
};
#endif // SOPLEX_WITH_CPPMPF

#endif // SOPLEX_WITH_BOOST

/**
 * @brief Templated tolerances class for arbitrary precision numeric types.
 *
 * This class stores numerical tolerances used throughout the SoPlex solver.
 * The template parameter R allows using different numeric types (double, long double,
 * or arbitrary precision types like Boost multiprecision).
 *
 * For high-precision solving, the class also supports storing epsilon values as
 * strings which can then be converted to any numeric type on access. This prevents
 * precision loss when storing values like 1e-600.
 *
 * For backward compatibility, a typedef `Tolerances = TolerancesBase<Real>` is provided.
 */
template <class R>
class TolerancesBase
{
private:

   //------------------------------------
   /**@name Data */
   ///@{
   /// default allowed additive zero: 1.0 + EPS_ZERO == 1.0
   R s_epsilon;
   /// epsilon for factorization
   R s_epsilon_factorization;
   /// epsilon for factorization update
   R s_epsilon_update;
   /// epsilon for pivot zero tolerance in factorization
   R s_epsilon_pivot;
   /// feasibility tolerance
   R s_feastol;
   /// optimality tolerance
   R s_opttol;
   /// floating point feasibility tolerance
   R s_floating_point_feastol;
   /// floating point optimality tolerance
   R s_floating_point_opttol;
   /// multiplier for fixed numbers that should change if s_epsilon changes
   R s_epsilon_multiplier;
   /// infinity threshold (R-typed for high precision)
   R s_infinity;

#ifdef SOPLEX_WITH_BOOST
   /// High-precision string storage for epsilon values (optional)
   /// When set via setEpsilonRational(), these strings store the exact value
   /// which can then be converted to any precision type on access
   std::string s_epsilon_str;
   std::string s_epsilon_factorization_str;
   std::string s_epsilon_update_str;
   std::string s_epsilon_pivot_str;
   std::string s_feastol_str;
   std::string s_opttol_str;
   bool s_rational_epsilon_set;           ///< true if rational epsilon strings are set
#endif
   ///@}

public:

   /// default constructor using precision-appropriate defaults
   explicit TolerancesBase()
      : s_epsilon(PrecisionTraits<R>::defaultEpsilon())
      , s_epsilon_factorization(PrecisionTraits<R>::defaultEpsilonFactorization())
      , s_epsilon_update(PrecisionTraits<R>::defaultEpsilonUpdate())
      , s_epsilon_pivot(PrecisionTraits<R>::defaultEpsilonPivot())
      , s_feastol(PrecisionTraits<R>::defaultFeastol())
      , s_opttol(PrecisionTraits<R>::defaultFeastol())
      , s_floating_point_feastol(PrecisionTraits<R>::defaultFeastol())
      , s_floating_point_opttol(PrecisionTraits<R>::defaultFeastol())
      , s_epsilon_multiplier(R(1))
      , s_infinity(PrecisionTraits<R>::defaultInfinity())
#ifdef SOPLEX_WITH_BOOST
      , s_epsilon_str("")
      , s_epsilon_factorization_str("")
      , s_epsilon_update_str("")
      , s_epsilon_pivot_str("")
      , s_feastol_str("")
      , s_opttol_str("")
      , s_rational_epsilon_set(false)
#endif
   {}

   //------------------------------------
   /**@name Access / modification */
   ///@{

   /// global zero epsilon
   R epsilon() const
   {
#ifdef SOPLEX_WITH_BOOST

      // Return string-converted value if high-precision epsilon is set
      if(s_rational_epsilon_set && !s_epsilon_str.empty())
         return StringToNumber<R>::convert(s_epsilon_str);

#endif
      return s_epsilon;
   }

   /// get epsilon in arbitrary precision type T (for high-precision operations)
   template <class T>
   T epsilonAs() const
   {
#ifdef SOPLEX_WITH_BOOST

      if(s_rational_epsilon_set && !s_epsilon_str.empty())
         return StringToNumber<T>::convert(s_epsilon_str);

#endif
      return T(s_epsilon);
   }

   /// set global zero epsilon
   void setEpsilon(const R& eps)
   {
      s_epsilon = eps;
      // Update epsilon multiplier for scaling
      R defaultEps = PrecisionTraits<R>::defaultEpsilon();

      if(defaultEps > R(0))
      {
         // Use sqrt of ratio to scale other tolerances proportionally
         R ratio = s_epsilon / defaultEps;

         if(ratio > R(0))
            s_epsilon_multiplier = SafeSqrt<R>::compute(ratio);
         else
            s_epsilon_multiplier = R(1);
      }
      else
      {
         s_epsilon_multiplier = R(1);
      }
   }

#ifdef SOPLEX_WITH_BOOST
   /// set epsilon from string value (for high-precision solving)
   /// This stores the exact string which is then converted to any precision on access
   void setEpsilonRational(const std::string& eps)
   {
      s_epsilon_str = eps;
      s_rational_epsilon_set = true;
   }

   /// set all epsilon values from strings (for high-precision solving)
   void setEpsilonsRational(const std::string& eps, const std::string& epsFactor,
                            const std::string& epsUpdate, const std::string& epsPivot)
   {
      s_epsilon_str = eps;
      s_epsilon_factorization_str = epsFactor;
      s_epsilon_update_str = epsUpdate;
      s_epsilon_pivot_str = epsPivot;
      s_rational_epsilon_set = true;
   }

   /// set epsilon factorization from string value (for high-precision solving)
   void setEpsilonFactorizationRational(const std::string& eps)
   {
      s_epsilon_factorization_str = eps;
      s_rational_epsilon_set = true;
   }

   /// set epsilon update from string value (for high-precision solving)
   void setEpsilonUpdateRational(const std::string& eps)
   {
      s_epsilon_update_str = eps;
      s_rational_epsilon_set = true;
   }

   /// set epsilon pivot from string value (for high-precision solving)
   void setEpsilonPivotRational(const std::string& eps)
   {
      s_epsilon_pivot_str = eps;
      s_rational_epsilon_set = true;
   }

   /// check if rational epsilons are set
   bool hasRationalEpsilons() const
   {
      return s_rational_epsilon_set;
   }

   /// get the stored epsilon string (for inspection)
   const std::string& epsilonRational() const
   {
      return s_epsilon_str;
   }

   /// get the stored factorization epsilon string
   const std::string& epsilonFactorizationRational() const
   {
      return s_epsilon_factorization_str;
   }

   /// get the stored update epsilon string
   const std::string& epsilonUpdateRational() const
   {
      return s_epsilon_update_str;
   }

   /// get the stored pivot epsilon string
   const std::string& epsilonPivotRational() const
   {
      return s_epsilon_pivot_str;
   }

   /// set feastol from string value (for high-precision solving)
   void setFeastolRational(const std::string& ftol)
   {
      s_feastol_str = ftol;
      s_rational_epsilon_set = true;
   }

   /// set opttol from string value (for high-precision solving)
   void setOpttolRational(const std::string& otol)
   {
      s_opttol_str = otol;
      s_rational_epsilon_set = true;
   }

   /// get the stored feastol string
   const std::string& feastolRational() const
   {
      return s_feastol_str;
   }

   /// get the stored opttol string
   const std::string& opttolRational() const
   {
      return s_opttol_str;
   }
#endif

   /// zero epsilon used in factorization
   R epsilonFactorization() const
   {
#ifdef SOPLEX_WITH_BOOST

      if(s_rational_epsilon_set && !s_epsilon_factorization_str.empty())
         return StringToNumber<R>::convert(s_epsilon_factorization_str);

#endif
      return s_epsilon_factorization;
   }

   /// set zero epsilon used in factorization
   void setEpsilonFactorization(const R& eps)
   {
      s_epsilon_factorization = eps;
   }

   /// zero epsilon used in factorization update
   R epsilonUpdate() const
   {
#ifdef SOPLEX_WITH_BOOST

      if(s_rational_epsilon_set && !s_epsilon_update_str.empty())
         return StringToNumber<R>::convert(s_epsilon_update_str);

#endif
      return s_epsilon_update;
   }

   /// set zero epsilon used in factorization update
   void setEpsilonUpdate(const R& eps)
   {
      s_epsilon_update = eps;
   }

   /// zero epsilon used in pivot
   R epsilonPivot() const
   {
#ifdef SOPLEX_WITH_BOOST

      if(s_rational_epsilon_set && !s_epsilon_pivot_str.empty())
         return StringToNumber<R>::convert(s_epsilon_pivot_str);

#endif
      return s_epsilon_pivot;
   }

   /// set zero epsilon used in pivot
   void setEpsilonPivot(const R& eps)
   {
      s_epsilon_pivot = eps;
   }

   /// global feasibility tolerance
   R feastol() const
   {
#ifdef SOPLEX_WITH_BOOST

      if(s_rational_epsilon_set && !s_feastol_str.empty())
         return StringToNumber<R>::convert(s_feastol_str);

#endif
      return s_feastol;
   }

   /// set global feasibility tolerance
   void setFeastol(const R& ftol)
   {
      s_feastol = ftol;
   }

   /// global optimality tolerance
   R opttol() const
   {
#ifdef SOPLEX_WITH_BOOST

      if(s_rational_epsilon_set && !s_opttol_str.empty())
         return StringToNumber<R>::convert(s_opttol_str);

#endif
      return s_opttol;
   }

   /// set global optimality tolerance
   void setOpttol(const R& otol)
   {
      s_opttol = otol;
   }

   /// floating point feasibility tolerance used within the solver
   R floatingPointFeastol() const
   {
      return s_floating_point_feastol;
   }

   /// set floating point feasibility tolerance used within the solver
   void setFloatingPointFeastol(const R& ftol)
   {
      s_floating_point_feastol = ftol;
   }

   /// floating point optimality tolerance used within the solver
   R floatingPointOpttol() const
   {
      return s_floating_point_opttol;
   }

   /// set floating point optimality tolerance used within the solver
   void setFloatingPointOpttol(const R& otol)
   {
      s_floating_point_opttol = otol;
   }

   /// Get marker value for sparse vector operations
   R marker() const
   {
      return std::numeric_limits<R>::min();
   }

   /// R-typed infinity threshold (replaces double-precision REALPARAM::INFTY)
   R infinity() const
   {
      return s_infinity;
   }

   /// set R-typed infinity threshold
   void setInfinity(const R& inf)
   {
      s_infinity = inf;
   }

   /// scale a value such that it remains unchanged at default epsilon, but is scaled with smaller epsilon values
   /// this is updated in setEpsilon()
   R scaleAccordingToEpsilon(const R& a) const
   {
      return s_epsilon_multiplier == R(1) ? a : a * s_epsilon_multiplier;
   }
   ///@}
};

/// Backward-compatible typedef: Tolerances uses Real (double/long double/float)
typedef TolerancesBase<Real> Tolerances;

// A generic version of spxAbs. It would be nice if we could replace spxAbs
// with std::abs. Currently there are different versions of spxAbs under
// compile time #if. It's better to make this an overloaded function. Even
// better, replace it by std::abs since types from boost/multiprecision would
// need no extra modification.
template <class R>
R spxAbs(R a)
{
   return abs(a);
}

// cmath means proper long double function gets called, e.g. for fabs -> fabsl.
// Documentation unclear for nextafterl, so the ifdef remains for that case.
#ifdef WITH_LONG_DOUBLE
// returns the next representable value after x in the direction of y
inline Real spxNextafter(Real x, Real y)
{
   return nextafterl(x, y);
}
#else
// returns the next representable value after x in the direction of y
inline Real spxNextafter(Real x, Real y)
{
#ifndef _MSC_VER
   return nextafter(x, y);
#else
   return _nextafter(x, y);
#endif
}
#endif

/// returns |a|
template <>
inline Real spxAbs(Real a)
{
   return fabs(a);
}

/// returns square root
inline Real spxSqrt(Real a)
{
   return std::sqrt(a);
}

/// returns max(|a|,|b|)
inline Real maxAbs(Real a, Real b)
{
   const Real absa = spxAbs(a);
   const Real absb = spxAbs(b);

   return absa > absb ? absa : absb;
}

/// returns (a-b) / max(|a|,|b|,1.0)
inline Real relDiff(Real a, Real b)
{
   return (a - b) / (maxAbs(a, b) > 1.0 ? maxAbs(a, b) : 1.0);
}

/// safe version of snprintf
inline int spxSnprintf(
   char*                 t,                  /**< target string */
   size_t                len,                /**< length of the string to copy */
   const char*           s,                  /**< source string */
   ...                                       /**< further parameters */
)
{
   va_list ap;
   int n;

   assert(t != nullptr);
   assert(len > 0);

   va_start(ap, s); /*lint !e826*/

#if defined(_WIN32) || defined(_WIN64)
   n = _vsnprintf(t, len, s, ap);
#else
   n = vsnprintf(t, len, s, ap); /*lint !e571*/
#endif
   va_end(ap);

   if(n < 0 || (size_t) n >= len)
   {
#ifndef NDEBUG

      if(n < 0)
      {
         SPX_MSG_ERROR(std::cerr << "vsnprintf returned " << n << " while reading: " << s << std::endl;)
      }

#endif
      t[len - 1] = '\0';
      n = (int) len - 1;
   }

   return n;
}

#ifdef SOPLEX_WITH_BOOST

using namespace boost::multiprecision;

#ifdef SOPLEX_WITH_GMP
template<boost::multiprecision::expression_template_option eto>
inline number<gmp_rational, eto> ldexp(number<gmp_rational, eto>, int exp)
{
   assert(false);
   return number<gmp_rational>();
}

template<boost::multiprecision::expression_template_option eto>
inline number<gmp_rational, eto> frexp(number<gmp_rational, eto>, int* exp)
{
   assert(false);
   return number<gmp_rational>();
}
#else
inline cpp_rational ldexp(cpp_rational r, int exp)
{
   assert(false);
   return cpp_rational();
}

inline cpp_rational frexp(cpp_rational, int* exp)
{
   assert(false);
   return cpp_rational();
}
#endif

// wrapped frexp function
template <typename T, boost::multiprecision::expression_template_option eto>
boost::multiprecision::number<T, eto> spxFrexp(boost::multiprecision::number<T, eto> y, int* exp)
{
   return frexp(y, exp);
}

// Overloaded spxLdexp
template <typename T, boost::multiprecision::expression_template_option eto>
boost::multiprecision::number<T> spxLdexp(boost::multiprecision::number<T, eto> x, int exp)
{
   return ldexp(x, exp);
}

// Overloaded function to return the square-root
template <typename T, expression_template_option ep>
number<T, ep> spxSqrt(number<T, ep> a)
{
   return sqrt(a);
}

// the nextafter function
template <typename T, expression_template_option eto>
number<T, eto> spxNextafter(number<T, eto> x,
                            number<T, eto> y)
{
   // Turns out that nextafter is not supported in the mpfr library? The mpfr
   // library does a different function named nextabove. Probably a
   // replacement? I've made an issue about this.
   // return nextafter(x,y);

   // @todo Temporarily, I'm returning 0
   assert(false);
   return 0;
}

// Returns the square root
template <typename T>
number<T> spxSqrt(number<T> a)
{
   return sqrt(a);
}

/// returns max(|a|,|b|)
template <typename T, expression_template_option et>
inline number<T, et> maxAbs(
   number<T, et> a, number<T, et> b)
{
   const auto absa = spxAbs(a);
   const auto absb = spxAbs(b);

   return absa > absb ? absa : absb;
}

template <typename T, expression_template_option et>
inline number<T, et> relDiff(number<T, et> a,
                             number<T, et> b)
{
   return (a - b) / (maxAbs(a, b) > 1.0 ? maxAbs(a, b) : 1.0);
}
#endif
using namespace soplex;

} // namespace soplex

// For the templated functions
#include "spxdefines.hpp"

#endif // _SPXDEFINES_H_
