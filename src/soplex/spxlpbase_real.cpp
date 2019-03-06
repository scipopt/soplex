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

/**@file  spxlpbase_real.hpp
 * @brief Saving LPs with R values in a form suitable for SoPlex.
 */

#include <assert.h>
#include <stdio.h>
#include <ctype.h>
#include <iostream>

#include "soplex/spxdefines.h"
#include "soplex/spxlpbase.h"
#include "soplex/spxout.h"
#include "soplex/mpsinput.h"
#include "soplex/exceptions.h"
#include "soplex/spxscaler.h"

// @todo #if, else
#include "boost/multiprecision/number.hpp"
#include <boost/multiprecision/mpfr.hpp>

namespace soplex{

using mpfr_float_50_eto = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<50>, boost::multiprecision::et_off>;

// @todo write boost
/// Changes objective vector to \p newObj.
template <>
void SPxLPBase<Real>::changeMaxObj(const VectorBase<Real>& newObj, bool scale)
{
  assert(maxObj().dim() == newObj.dim());
  if( scale )
    {
      assert(_isScaled);
      assert(lp_scaler);
      LPColSetBase<Real>::maxObj_w().scaleAssign(LPColSetBase<Real>::scaleExp.get_const_ptr(), newObj);
    }
  else
    LPColSetBase<Real>::maxObj_w() = newObj;
  assert(isConsistent());
}

  // @todo: uncomment this
// template <>
// void SPxLPBase<mpfr_float_50_eto>::changeMaxObj(const VectorBase<mpfr_float_50_eto>& newObj, bool scale)
// {
//   assert(maxObj().dim() == newObj.dim());
//   if( scale )
//     {
//       assert(_isScaled);
//       assert(lp_scaler);
//       LPColSetBase<mpfr_float_50_eto>::maxObj_w().scaleAssign(LPColSetBase<mpfr_float_50_eto>::scaleExp.get_const_ptr(), newObj);
//     }
//   else
//     LPColSetBase<mpfr_float_50_eto>::maxObj_w() = newObj;
//   assert(isConsistent());
// }

/// Changes vector of lower bounds to \p newLower.
template <>
void SPxLPBase<Real>::changeLower(const VectorBase<Real>& newLower, bool scale)
{
  assert(lower().dim() == newLower.dim());
  if( scale )
    {
      assert(_isScaled);
      assert(lp_scaler);
      LPColSetBase<Real>::lower_w().scaleAssign(LPColSetBase<Real>::scaleExp.get_const_ptr(), newLower, true);
    }
  else
    LPColSetBase<Real>::lower_w() = newLower;
  assert(isConsistent());
}

// template <>
// void SPxLPBase<mpfr_float_50_eto>::changeLower(const VectorBase<mpfr_float_50_eto>& newLower, bool scale)
// {
//   assert(lower().dim() == newLower.dim());
//   if( scale )
//     {
//       assert(_isScaled);
//       assert(lp_scaler);
//       LPColSetBase<mpfr_float_50_eto>::lower_w().scaleAssign(LPColSetBase<mpfr_float_50_eto>::scaleExp.get_const_ptr(), newLower, true);
//     }
//   else
//     LPColSetBase<mpfr_float_50_eto>::lower_w() = newLower;
//   assert(isConsistent());
// }


/// Changes vector of upper bounds to \p newUpper.
template <>
void SPxLPBase<Real>::changeUpper(const VectorBase<Real>& newUpper, bool scale)
{
  assert(upper().dim() == newUpper.dim());
  if( scale )
    {
      assert(_isScaled);
      assert(lp_scaler);
      LPColSetBase<Real>::upper_w().scaleAssign(LPColSetBase<Real>::scaleExp.get_const_ptr(), newUpper, true);
    }
  else
    LPColSetBase<Real>::upper_w() = newUpper;
  assert(isConsistent());
}


// /// Changes vector of upper bounds to \p newUpper.
// template <>
// void SPxLPBase<mpfr_float_50_eto>::changeUpper(const VectorBase<mpfr_float_50_eto>& newUpper, bool scale)
// {
//   assert(upper().dim() == newUpper.dim());
//   if( scale )
//     {
//       assert(_isScaled);
//       assert(lp_scaler);
//       LPColSetBase<mpfr_float_50_eto>::upper_w().scaleAssign(LPColSetBase<mpfr_float_50_eto>::scaleExp.get_const_ptr(), newUpper, true);
//     }
//   else
//     LPColSetBase<mpfr_float_50_eto>::upper_w() = newUpper;
//   assert(isConsistent());
// }

// @todo; the boost version of the following function
/// Changes left hand side vector for constraints to \p newLhs.
template <>
void SPxLPBase<Real>::changeLhs(const VectorBase<Real>& newLhs, bool scale)
{
  assert(lhs().dim() == newLhs.dim());
  if( scale )
    {
      assert(_isScaled);
      assert(lp_scaler);
      LPRowSetBase<Real>::lhs_w().scaleAssign(LPRowSetBase<Real>::scaleExp.get_const_ptr(), newLhs);
    }
  else
    LPRowSetBase<Real>::lhs_w() = newLhs;
  assert(isConsistent());
}

// @todo: the boost version of the following function
/// Changes right hand side vector for constraints to \p newRhs.
template <>
void SPxLPBase<Real>::changeRhs(const VectorBase<Real>& newRhs, bool scale)
{
  assert(rhs().dim() == newRhs.dim());
  if( scale )
    {
      assert(_isScaled);
      assert(lp_scaler);
      LPRowSetBase<Real>::rhs_w().scaleAssign(LPRowSetBase<Real>::scaleExp.get_const_ptr(), newRhs);
    }
  else
    LPRowSetBase<Real>::rhs_w() = newRhs;
  assert(isConsistent());
}

  // @todo: fix this file and remove unnecessary functions.
/// Changes LP element (\p i, \p j) to \p val. \p scale determines whether the new data should be scaled
template <>
void SPxLPBase<mpfr_float_50_eto>::changeElement(int i, int j, const mpfr_float_50_eto &val, bool scale)
{
   if (i < 0 || j < 0)
      return;

   SVectorBase<mpfr_float_50_eto> &row = rowVector_w(i);
   SVectorBase<mpfr_float_50_eto> &col = colVector_w(j);

   if( isNotZero(val) )
   {
     mpfr_float_50_eto newVal;

      if (scale)
      {
         assert(_isScaled);
         assert(lp_scaler);
         newVal = lp_scaler->scaleElement(*this, i, j, val);
      }
      else
         newVal = val;

      if (row.pos(j) >= 0 && col.pos(i) >= 0)
      {
         row.value(row.pos(j)) = newVal;
         col.value(col.pos(i)) = newVal;
      }
      else
      {
         LPRowSetBase<mpfr_float_50_eto>::add2(i, 1, &j, &newVal);
         LPColSetBase<mpfr_float_50_eto>::add2(j, 1, &i, &newVal);
      }
  }
   else if (row.pos(j) >= 0 && col.pos(i) >= 0)
   {
      row.remove(row.pos(j));
      col.remove(col.pos(i));
   }

   assert(isConsistent());
}

}
