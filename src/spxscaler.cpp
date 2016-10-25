/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  spxscaler.cpp
 * @brief LP scaling base class.
 */

#define BITSHIFTSCALING

#ifdef BITSHIFTSCALING
#include <cmath>
#endif

#include <iostream>
#include <assert.h>

#include "spxscaler.h"
#include "spxlp.h"
#include "dsvector.h"
#include "dvector.h"

namespace soplex
{

std::ostream& operator<<(std::ostream& s, const SPxScaler& sc)
{
   s << sc.getName() << " scaler:" << std::endl;
   s << "colscale = [ ";
   for(int ci = 0; ci < sc.m_colscaleExp.size(); ++ci )
      s << sc.m_colscaleExp[ci] << " ";
   s << "]" << std::endl;

   s << "rowscale = [ ";
   for(int ri = 0; ri < sc.m_rowscaleExp.size(); ++ri )
      s << sc.m_rowscaleExp[ri] << " ";
   s << "]" << std::endl;

   return s;
}

SPxScaler::SPxScaler(
   const char* name, 
   bool        colFirst, 
   bool        doBoth,
   SPxOut*     outstream)
   : m_name(name)
   , m_colFirst(colFirst)
   , m_doBoth(doBoth)
   , spxout(outstream)
{
   assert(SPxScaler::isConsistent());
}

SPxScaler::SPxScaler(const SPxScaler& old)
   : m_name(old.m_name)
   , m_colscaleExp(old.m_colscaleExp)
   , m_rowscaleExp(old.m_rowscaleExp)
   , m_colFirst(old.m_colFirst)
   , m_doBoth(old.m_doBoth)
   , spxout(old.spxout)
{
   assert(SPxScaler::isConsistent());
}

SPxScaler::~SPxScaler()
{
   m_name = 0;
}

SPxScaler& SPxScaler::operator=(const SPxScaler& rhs)
{
   if (this != &rhs)
   {
      m_name     = rhs.m_name;
      m_colscaleExp = rhs.m_colscaleExp;
      m_rowscaleExp = rhs.m_rowscaleExp;
      m_colFirst = rhs.m_colFirst;
      m_doBoth   = rhs.m_doBoth;
      spxout     = rhs.spxout;

      assert(SPxScaler::isConsistent());
   }
   return *this;
}


const char* SPxScaler::getName() const
{

   return m_name;
}

void SPxScaler::setOrder(bool colFirst)
{

   m_colFirst = colFirst;
}

void SPxScaler::setBoth(bool both)
{

   m_doBoth = both;
}

void SPxScaler::setRealParam(Real param, const char* name)
{}

void SPxScaler::setIntParam(int param, const char* name)
{}

void SPxScaler::setup(SPxLP& lp)
{

   assert(lp.isConsistent());

   m_colscaleExp.reSize(lp.nCols());
   m_rowscaleExp.reSize(lp.nRows());

   int i;

   for(i = 0; i < lp.nCols(); ++i )
      m_colscaleExp[i] = 0.0;

   for(i = 0; i < lp.nRows(); ++i )
      m_rowscaleExp[i] = 0.0;
}

#if 0
/** This function is used by computeScaleVecs and has to be overridden.
 */
Real SPxScaler::computeScale(Real /*mini*/, Real /*maxi*/) const
{

   return 1.0;
}

Real SPxScaler::computeScalingVecs(
   const SVSet*           vecset, 
   const DataArray<int>& coScaleExp,
   DataArray<int>&       scaleExp)
{

   Real pmax = 0.0;

   for( int i = 0; i < vecset->num(); ++i )
   {
      const SVector& vec = (*vecset)[i];

      Real maxi = 0.0;
      Real mini = infinity;

      for( int j = 0; j < vec.size(); ++j )
      {
         Real x = spxAbs(vec.value(j) * spxLdexp(1.0, coScaleExp[vec.index(j)]));

         if (!isZero(x))
         {
            if (x > maxi)
               maxi = x;
            if (x < mini)
               mini = x;
         }
      }
      // empty rows/cols are possible
      if (mini == infinity || maxi == 0.0)
      {
         mini = 1.0;
         maxi = 1.0;
      }
      assert(mini < infinity);
      assert(maxi > 0.0);

      //scaleval[i] = 1.0 / computeScale(mini, maxi);
      frexp(1.0 / computeScale(mini, maxi), &(scaleExp[i]));
      scaleExp[i] -= 1;

      Real p = maxi / mini;

      if (p > pmax)
         pmax = p;
   }
   return pmax;
}
#endif
void SPxScaler::applyScaling(SPxLP& lp)
{
   for( int i = 0; i < lp.nRows(); ++i )
   {
      SVector& vec = lp.rowVector_w(i);
#ifdef BITSHIFTSCALING
      int exp1;
      //frexp(m_rowscaleExp[i], &exp2);
      int exp2 = m_rowscaleExp[i];

      for( int j = 0; j < vec.size(); ++j)
      {
         //frexp(m_colscaleExp[vec.index(j)], &exp1);
         exp1 = m_colscaleExp[vec.index(j)];
         vec.value(j) = spxLdexp(vec.value(j), exp1 + exp2);
      }
      if (lp.rhs(i) < infinity)
      {
         lp.rhs_w(i) = spxLdexp(lp.rhs_w(i), exp2);
      }
      if (lp.lhs(i) > -infinity)
      {
         lp.lhs_w(i) = spxLdexp(lp.lhs_w(i), exp2);

      }
#else
      for( int j = 0; j < vec.size(); ++j)
         vec.value(j) *= m_colscaleExp[vec.index(j)] * m_rowscaleExp[i];

      if (lp.rhs(i) < infinity)
         lp.rhs_w(i) *= m_rowscaleExp[i];
      if (lp.lhs(i) > -infinity)
         lp.lhs_w(i) *= m_rowscaleExp[i];
#endif
   }

   for( int i = 0; i < lp.nCols(); ++i )
   {
      SVector& vec = lp.colVector_w(i);
#ifdef BITSHIFTSCALING
      int exp1;
      int exp2 = m_colscaleExp[i];

      for( int j = 0; j < vec.size(); ++j)
      {
         exp1 = m_rowscaleExp[vec.index(j)];
         vec.value(j) = spxLdexp(vec.value(j), exp1 + exp2);
      }

      lp.maxObj_w(i) = spxLdexp(lp.maxObj_w(i), exp2);

      if (lp.upper(i) < infinity)
      {
         lp.upper_w(i) = spxLdexp(lp.upper_w(i), -exp2);
      }
      if (lp.lower(i) > -infinity)
      {
         lp.lower_w(i) = spxLdexp(lp.lower_w(i), -exp2);
      }
#else
      for( int j = 0; j < vec.size(); ++j)
         vec.value(j) *= m_rowscaleExp[vec.index(j)] * m_colscaleExp[i];

      lp.maxObj_w(i) *= m_colscaleExp[i];

      if (lp.upper(i) < infinity)
         lp.upper_w(i) /= m_colscaleExp[i];
      if (lp.lower(i) > -infinity)
         lp.lower_w(i) /= m_colscaleExp[i];
#endif
   }
   assert(lp.isConsistent());
}

/// unscale SPxLP
void SPxScaler::unscale(SPxLPBase<Real>& lp)
{
   for( int i = 0; i < lp.nRows(); ++i )
   {
      SVector& vec = lp.rowVector_w(i);

      int exp1;
      int exp2 = m_rowscaleExp[i];

      for( int j = 0; j < vec.size(); ++j)
      {
         exp1 = m_colscaleExp[vec.index(j)];
         vec.value(j) = spxLdexp(vec.value(j), -exp1 - exp2);
      }
      if (lp.rhs(i) < infinity)
      {
         lp.rhs_w(i) = spxLdexp(lp.rhs_w(i), -exp2);
      }
      if (lp.lhs(i) > -infinity)
      {
         lp.lhs_w(i) = spxLdexp(lp.lhs_w(i), -exp2);
      }
   }
   for( int i = 0; i < lp.nCols(); ++i )
   {
      SVector& vec = lp.colVector_w(i);

      int exp1;
      int exp2 = m_colscaleExp[i];

      for( int j = 0; j < vec.size(); ++j)
      {
         exp1 = m_rowscaleExp[vec.index(j)];
         vec.value(j) = spxLdexp(vec.value(j), -exp1 - exp2);
      }

      lp.maxObj_w(i) = spxLdexp(lp.maxObj_w(i), -exp2);

      if (lp.upper(i) < infinity)
      {
         lp.upper_w(i) = spxLdexp(lp.upper_w(i), exp2);
      }
      if (lp.lower(i) > -infinity)
      {
         lp.lower_w(i) = spxLdexp(lp.lower_w(i), exp2);
      }
   }
   assert(lp.isConsistent());
}

/// returns scaling factor for column \p i
Real SPxScaler::getColScaleExp(int i)
{
   return m_colscaleExp[i];
}

/// returns scaling factor for row \p i
Real SPxScaler::getRowScaleExp(int i)
{
   return m_rowscaleExp[i];
}


/// Gets unscaled column \p i
void SPxScaler::getColUnscaled(const SPxLP& lp, int i, SVector& vec) const
{
   assert(i < lp.nCols());
   assert(i >= 0);

   vec = lp.LPColSet::colVector(i);

   int exp1;
   int exp2 = m_colscaleExp[i];

   for( int j = 0; j < vec.size(); j++ )
   {
      exp1 = m_rowscaleExp[vec.index(j)];
      vec.value(j) = spxLdexp(vec.value(j), -exp1 - exp2);
   }
}


/// returns unscaled upper bound \p i
Real SPxScaler::upperUnscaled(const SPxLPBase<Real>& lp, int i) const
{
   assert(i < lp.nCols());
   assert(i >= 0);

   int exp = m_colscaleExp[i];

   if( lp.LPColSet::upper(i) < infinity )
   {
      return spxLdexp(lp.LPColSet::upper(i) , exp - 1);
   }
   else
   {
      return lp.LPColSet::upper(i);
   }
}


/// gets unscaled upper bound vector
void SPxScaler::getUpperUnscaled(const SPxLPBase<Real>& lp, DVector& vec) const
{
   int exp;

   for( int i = 0; i < lp.LPColSet::upper().dim(); i++)
   {
      exp = m_colscaleExp[i];
      vec[i] = spxLdexp(lp.LPColSet::upper()[i], exp - 1);
   }
}


/// returns unscaled upper bound vector of LP \lp
Real SPxScaler::lowerUnscaled(const SPxLPBase<Real>& lp, int i) const
{
   assert(i < lp.nCols());
   assert(i >= 0);

   int exp = m_colscaleExp[i];

   if( lp.LPColSet::lower(i) > -infinity )
   {
      return spxLdexp(lp.LPColSet::lower(i), exp - 1);
   }
   else
   {
      return lp.LPColSet::lower(i);
   }
}


/// returns unscaled lower bound vector of LP \lp
void SPxScaler::getLowerUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const
{
   int exp;

   for( int i = 0; i < lp.LPColSet::lower().dim(); i++)
   {
      exp = m_colscaleExp[i];
      vec[i] = spxLdexp(lp.LPColSet::lower()[i], exp - 1);
   }
}

/// returns unscaled objective function coefficient of \p i
Real SPxScaler::maxObjUnscaled(const SPxLPBase<Real>& lp, int i) const
{
   assert(i < lp.nCols());
   assert(i >= 0);

   int exp = m_colscaleExp[i];

   return spxLdexp(lp.LPColSet::maxObj(i) , -exp + 1);
}


/// gets unscaled objective function coefficient of \p i
void SPxScaler::getMaxObjUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const
{
   int exp;

   for( int i = 0; i < lp.LPColSet::maxObj().dim(); i++)
   {
      exp = m_colscaleExp[i];
      vec[i] = spxLdexp(lp.LPColSet::maxObj()[i], -exp + 1);
   }
}

/// returns unscaled row \p i
void SPxScaler::getRowUnscaled(const SPxLP& lp, int i, SVector& vec) const
{
   assert(i < lp.nRows());
   assert(i >= 0);

   int exp1;
   int exp2 = m_rowscaleExp[i];

   for( int j = 0; j < vec.size(); j++ )
   {
      exp1 = m_colscaleExp[vec.index(j)];
      vec.value(j) = spxLdexp(vec.value(j), -exp1 - exp2);
   }
}

/// returns unscaled right hand side \p i
Real SPxScaler::rhsUnscaled(const SPxLPBase<Real>& lp, int i) const
{
   assert(i < lp.nRows());
   assert(i >= 0);

   int exp = m_rowscaleExp[i];

   if( lp.LPRowSet::rhs(i) < infinity )
   {
      return spxLdexp(lp.LPRowSet::rhs(i) , -exp + 1);
   }
   else
   {
      return lp.LPRowSet::rhs(i);
   }
}


/// gets unscaled right hand side vector
void SPxScaler::getRhsUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const
{
   int exp;

   for( int i = 0; i < lp.LPRowSet::rhs().dim(); i++)
   {
      exp = m_rowscaleExp[i];
      vec[i] = spxLdexp(lp.LPRowSet::rhs()[i], -exp + 1);
   }
}


/// returns unscaled left hand side \p i of LP \lp
Real SPxScaler::lhsUnscaled(const SPxLPBase<Real>& lp, int i) const
{
   assert(i < lp.nRows());
   assert(i >= 0);

   int exp = m_rowscaleExp[i];

   if( lp.LPRowSet::lhs(i) > -infinity )
   {
      return spxLdexp(lp.LPRowSet::lhs(i) , -exp + 1);
   }
   else
   {
      return lp.LPRowSet::lhs(i);
   }
}

/// returns unscaled left hand side vector of LP \lp
void SPxScaler::getLhsUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const
{
   int exp;

   for( int i = 0; i < lp.LPRowSet::lhs().dim(); i++)
   {
      exp = m_rowscaleExp[i];
      vec[i] = spxLdexp(lp.LPRowSet::lhs()[i], -exp + 1);
   }
}

void SPxScaler::unscalePrimal(Vector& x) const
{

   assert(x.dim() == m_colscaleExp.size());
#ifdef BITSHIFTSCALING
   int exp1;
   for( int j = 0; j < x.dim(); ++j )
   {
      exp1 = m_colscaleExp[j];
      x[j] = spxLdexp(x[j], exp1);
   }
#else
   for(int j = 0; j < x.dim(); ++j)
      x[j] *= m_colscaleExp[j];
#endif
}

void SPxScaler::unscaleSlacks(Vector& s) const
{

   assert(s.dim() == m_rowscaleExp.size());
#ifdef BITSHIFTSCALING
   int exp1;
   for( int i = 0; i < s.dim(); ++i )
   {
      exp1 = m_rowscaleExp[i];
      s[i] = spxLdexp(s[i], -exp1);
   }
#else
   for(int i = 0; i < s.dim(); ++i)
      s[i] /= m_rowscaleExp[i];
#endif
}

void SPxScaler::unscaleDual(Vector& pi) const
{

   assert(pi.dim() == m_rowscaleExp.size());
#ifdef BITSHIFTSCALING
   int exp1;
   for( int i = 0; i < pi.dim(); ++i )
   {
      exp1 = m_rowscaleExp[i];
      pi[i] = spxLdexp(pi[i], exp1);
   }
#else
   for(int i = 0; i < pi.dim(); ++i)
      pi[i] *= m_rowscaleExp[i];
#endif
}

void SPxScaler::unscaleRedCost(Vector& r) const
{

   assert(r.dim() == m_colscaleExp.size());
#ifdef BITSHIFTSCALING
   int exp1;
   for( int j = 0; j < r.dim(); ++j )
   {
      exp1 = m_colscaleExp[j];
      r[j] = spxLdexp(r[j], -exp1);
   }
#else
   for(int j = 0; j < r.dim(); ++j)
      r[j] /= m_colscaleExp[j];
#endif
}

Real SPxScaler::minAbsColscale() const
{

   Real mini = infinity;

   for( int i = 0; i < m_colscaleExp.size(); ++i )
      if( spxAbs(spxLdexp(1.0, m_colscaleExp[i])) < mini )
         mini = spxAbs(spxLdexp(1.0, m_colscaleExp[i]));

   return mini;
}

Real SPxScaler::maxAbsColscale() const
{

   Real maxi = 0.0;

   for( int i = 0; i < m_colscaleExp.size(); ++i )
      if( spxAbs(spxLdexp(1.0, m_colscaleExp[i])) > maxi )
         maxi = spxAbs(spxLdexp(1.0, m_colscaleExp[i]));


   return maxi;
}

Real SPxScaler::minAbsRowscale() const
{

   Real mini = infinity;

   for( int i = 0; i < m_rowscaleExp.size(); ++i )
      if( spxAbs(spxLdexp(1.0, m_rowscaleExp[i])) < mini )
         mini = spxAbs(spxLdexp(1.0, m_rowscaleExp[i]));

   return mini;
}

Real SPxScaler::maxAbsRowscale() const
{

   Real maxi = 0.0;

   for( int i = 0; i < m_rowscaleExp.size(); ++i )
      if( spxAbs(spxLdexp(1.0, m_rowscaleExp[i])) > maxi )
         maxi = spxAbs(spxLdexp(1.0, m_rowscaleExp[i]));

   return maxi;
}

/** \f$\max_{j\in\mbox{ cols}}
 *   \left(\frac{\max_{i\in\mbox{ rows}}|a_ij|}
 *              {\min_{i\in\mbox{ rows}}|a_ij|}\right)\f$
 */
Real SPxScaler::maxColRatio(const SPxLP& lp) const
{

   Real pmax = 0.0;

   for(int i = 0; i < lp.nCols(); ++i )
   {
      const SVector& vec  = lp.colVector(i);
      Real           mini = infinity;
      Real           maxi = 0.0;

      for(int j = 0; j < vec.size(); ++j)
      {
         Real x = spxAbs(vec.value(j));

         if (x < mini)
            mini = x;
         if (x > maxi)
            maxi = x;
      }
      Real p = maxi / mini;

      if (p > pmax)
         pmax = p;
   }
   return pmax;
}

/** \f$\max_{i\in\mbox{ rows}}
 *   \left(\frac{\max_{j\in\mbox{ cols}}|a_ij|}
 *              {\min_{j\in\mbox{ cols}}|a_ij|}\right)\f$
 */
Real SPxScaler::maxRowRatio(const SPxLP& lp) const
{

   Real pmax = 0.0;

   for(int i = 0; i < lp.nRows(); ++i )
   {
      const SVector& vec  = lp.rowVector(i);
      Real           mini = infinity;
      Real           maxi = 0.0;

      for(int j = 0; j < vec.size(); ++j)
      {
         Real x = spxAbs(vec.value(j));

         if (x < mini)
            mini = x;
         if (x > maxi)
            maxi = x;
      }
      Real p = maxi / mini;

      if (p > pmax)
         pmax = p;
   }
   return pmax;
}

bool SPxScaler::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS

   return m_colscaleExpExp.isConsistent() && m_rowscaleExpExp.isConsistent();
#else
   return true;
#endif
}

} // namespace soplex
