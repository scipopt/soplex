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

#ifndef _SOPLEX_STABLE_SUM_H_
#define _SOPLEX_STABLE_SUM_H_

#include <type_traits>

namespace soplex {

template <typename REAL>
class StableSum
{
   typename std::remove_const<REAL>::type sum = 0;

 public:
   StableSum() = default;
   StableSum( const REAL& init ) : sum( init ) {}

   void operator+=( const REAL& input )
   {
      sum += input;
   }

   void operator-=( const REAL& input )
   {
      sum -= input;
   }

   operator typename std::remove_const<REAL>::type() const
   {
      return sum;
   }
};

template <>
class StableSum<double>
{
   double sum = 0;
   double c = 0;

 public:
   StableSum() = default;
   StableSum( double init ) : sum( init ), c( 0 ) {}

   void operator+=( double input )
   {
#if _MSC_VER || __INTEL_COMPILER
      #pragma float_control( precise, on )
#endif

      double t = sum + input;
      double z = t - sum;
      double y = ( sum - ( t - z ) ) + ( input - z );
      c += y;

      sum = t;
   }

   void operator-=( double input )
   {
      (*this) += -input;
   }

   operator double() const
   {
      return sum + c;
   }
};
}

#endif
