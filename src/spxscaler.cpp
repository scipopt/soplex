/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2015 Konrad-Zuse-Zentrum                            */
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

namespace soplex
{

std::ostream& operator<<(std::ostream& s, const SPxScaler& sc)
{
   s << sc.getName() << " scaler:" << std::endl;
   s << "colscale = [ ";
   for(int ci = 0; ci < sc.m_colscale.size(); ++ci )
      s << sc.m_colscale[ci] << " ";
   s << "]" << std::endl;

   s << "rowscale = [ ";
   for(int ri = 0; ri < sc.m_rowscale.size(); ++ri )
      s << sc.m_rowscale[ri] << " ";
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
   , m_colscale(old.m_colscale)
   , m_rowscale(old.m_rowscale)
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
      m_colscale = rhs.m_colscale;
      m_rowscale = rhs.m_rowscale;
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

void SPxScaler::setup(SPxLP& lp)
{

   assert(lp.isConsistent());

   m_colscale.reSize(lp.nCols());
   m_rowscale.reSize(lp.nRows());

   int i;

   for(i = 0; i < lp.nCols(); ++i )
      m_colscale[i] = 1.0;

   for(i = 0; i < lp.nRows(); ++i )
      m_rowscale[i] = 1.0;
}

/** This function is used by computeScaleVecs and has to be overridden.
 */
Real SPxScaler::computeScale(Real /*mini*/, Real /*maxi*/) const
{

   return 1.0;
}

Real SPxScaler::computeScalingVecs(
   const SVSet*           vecset, 
   const DataArray<Real>& coScaleval, 
   DataArray<Real>&       scaleval) 
{

   Real pmax = 0.0;

   for(int i = 0; i < vecset->num(); ++i )
   {
      const SVector& vec = (*vecset)[i];

      Real maxi = 0.0;
      Real mini = infinity;

      for( int j = 0; j < vec.size(); ++j)
      {
         Real x = spxAbs(vec.value(j) * coScaleval[vec.index(j)]);

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

      scaleval[i] = 1.0 / computeScale(mini, maxi);

      Real p = maxi / mini;

      if (p > pmax)
         pmax = p;
   }
   return pmax;
}

void SPxScaler::applyScaling(SPxLP& lp)
{

   int i;

   for(i = 0; i < lp.nRows(); ++i )
   {
      SVector& vec = lp.rowVector_w(i);
#ifdef BITSHIFTSCALING
      int exp1,exp2;
      for( int j = 0; j < vec.size(); ++j)
      {
         spxFrexp(m_colscale[vec.index(j)], &exp1);
         spxFrexp(m_rowscale[i], &exp2);
         vec.value(j) = spxLdexp(vec.value(j), exp1 + exp2 - 2);
      }
      if (lp.rhs(i) < infinity)
      {
         spxFrexp(m_rowscale[i], &exp1);
         lp.rhs_w(i) = spxLdexp(lp.rhs_w(i), exp1 - 1);
      }
      if (lp.lhs(i) > -infinity)
      {
         spxFrexp(m_rowscale[i], &exp1);
         lp.lhs_w(i) = spxLdexp(lp.lhs_w(i), exp1 - 1);
      }
#else
      for( int j = 0; j < vec.size(); ++j)
         vec.value(j) *= m_colscale[vec.index(j)] * m_rowscale[i];

      if (lp.rhs(i) < infinity)
         lp.rhs_w(i) *= m_rowscale[i];
      if (lp.lhs(i) > -infinity)
         lp.lhs_w(i) *= m_rowscale[i];
#endif
   }
   for(i = 0; i < lp.nCols(); ++i )
   {
      SVector& vec = lp.colVector_w(i);
#ifdef BITSHIFTSCALING
      int exp1,exp2;
      for( int j = 0; j < vec.size(); ++j)
      {
         spxFrexp(m_rowscale[vec.index(j)], &exp1);
         spxFrexp(m_colscale[i], &exp2);
         vec.value(j) = spxLdexp(vec.value(j), exp1 + exp2 - 2);
      }

      spxFrexp(m_colscale[i], &exp1);
      lp.maxObj_w(i) = spxLdexp(lp.maxObj_w(i), exp1 - 1);

      if (lp.upper(i) < infinity)
      {
         spxFrexp(m_colscale[i], &exp1);
         lp.upper_w(i) = spxLdexp(lp.upper_w(i), -exp1 + 1);
      }
      if (lp.lower(i) > -infinity)
      {
         spxFrexp(m_colscale[i], &exp1);
         lp.lower_w(i) = spxLdexp(lp.lower_w(i), -exp1 + 1);
      }
#else
      for( int j = 0; j < vec.size(); ++j)
         vec.value(j) *= m_rowscale[vec.index(j)] * m_colscale[i];

      lp.maxObj_w(i) *= m_colscale[i];

      if (lp.upper(i) < infinity)
         lp.upper_w(i) /= m_colscale[i];
      if (lp.lower(i) > -infinity)
         lp.lower_w(i) /= m_colscale[i];
#endif
   }
   assert(lp.isConsistent());
}

void SPxScaler::unscalePrimal(Vector& x) const
{

   assert(x.dim() == m_colscale.size());
#ifdef BITSHIFTSCALING
   int exp1;
   for(int j = 0; j < x.dim(); ++j)
   {
      spxFrexp(m_colscale[j], &exp1);
      x[j] = spxLdexp(x[j], exp1 - 1);
   }
#else
   for(int j = 0; j < x.dim(); ++j)
      x[j] *= m_colscale[j];
#endif
}

void SPxScaler::unscaleSlacks(Vector& s) const
{

   assert(s.dim() == m_rowscale.size());
#ifdef BITSHIFTSCALING
   int exp1;
   for(int i = 0; i < s.dim(); ++i)
   {
      spxFrexp(m_rowscale[i], &exp1);
      s[i] = spxLdexp(s[i], -exp1 + 1);
   }
#else
   for(int i = 0; i < s.dim(); ++i)
      s[i] /= m_rowscale[i];
#endif
}

void SPxScaler::unscaleDual(Vector& pi) const
{

   assert(pi.dim() == m_rowscale.size());
#ifdef BITSHIFTSCALING
   int exp1;
   for(int i = 0; i < pi.dim(); ++i)
   {
      spxFrexp(m_rowscale[i], &exp1);
      pi[i] = spxLdexp(pi[i], exp1 - 1);
   }
#else
   for(int i = 0; i < pi.dim(); ++i)
      pi[i] *= m_rowscale[i];
#endif
}

void SPxScaler::unscaleRedCost(Vector& r) const
{

   assert(r.dim() == m_colscale.size());
#ifdef BITSHIFTSCALING
   int exp1;
   for(int j = 0; j < r.dim(); ++j)
   {
      spxFrexp(m_colscale[j], &exp1);
      r[j] = spxLdexp(r[j], -exp1 + 1);
   }
#else
   for(int j = 0; j < r.dim(); ++j)
      r[j] /= m_colscale[j];
#endif
}

Real SPxScaler::minAbsColscale() const
{

   Real mini = infinity;

   for(int i = 0; i < m_colscale.size(); ++i)
      if (spxAbs(m_colscale[i]) < mini)
         mini = spxAbs(m_colscale[i]);
#ifdef BITSHIFTSCALING
   int exp;
   spxFrexp(mini, &exp);
   mini = spxLdexp(2.0, exp - 1);
#endif
   return mini;
}

Real SPxScaler::maxAbsColscale() const
{

   Real maxi = 0.0;

   for(int i = 0; i < m_colscale.size(); ++i)
      if (spxAbs(m_colscale[i]) > maxi)
         maxi = spxAbs(m_colscale[i]);

#ifdef BITSHIFTSCALING
   int exp;
   spxFrexp(maxi, &exp);
   maxi = spxLdexp(2.0, exp - 1);
#endif
   return maxi;
}

Real SPxScaler::minAbsRowscale() const
{

   Real mini = infinity;

   for(int i = 0; i < m_rowscale.size(); ++i)
      if (spxAbs(m_rowscale[i]) < mini)
         mini = spxAbs(m_rowscale[i]);
#ifdef BITSHIFTSCALING
   int exp;
   spxFrexp(mini, &exp);
   mini = spxLdexp(2.0, exp - 1);
#endif
   return mini;
}

Real SPxScaler::maxAbsRowscale() const
{

   Real maxi = 0.0;

   for(int i = 0; i < m_rowscale.size(); ++i)
      if (spxAbs(m_rowscale[i]) > maxi)
         maxi = spxAbs(m_rowscale[i]);
#ifdef BITSHIFTSCALING
   int exp;
   spxFrexp(maxi, &exp);
   maxi = spxLdexp(2.0, exp - 1);
#endif
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

   return m_colscale.isConsistent() && m_rowscale.isConsistent();
#else
   return true;
#endif
}

} // namespace soplex
