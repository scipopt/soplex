/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxscaler.cpp,v 1.5 2003/01/05 19:03:17 bzfkocht Exp $"

/**@file  spxscaler.cpp
 * @brief LP scaling base class.
 */
#include <assert.h>

#include "spxscaler.h"
#include "spxlp.h"

namespace soplex
{

std::ostream& operator<<(std::ostream& s, const SPxScaler& sc)
{
   s << sc.getName() << " scaler:" << std::endl;
   s << "colscale = [ ";
   for(int i = 0; i < sc.m_colscale.size(); ++i )
      s << sc.m_colscale[i] << " ";
   s << "]" << std::endl;
      
   s << "rowscale = [ ";
   for(int i = 0; i < sc.m_rowscale.size(); ++i )
      s << sc.m_rowscale[i] << " ";
   s << "]" << std::endl;

   return s;
}

SPxScaler::SPxScaler(
   const char* name, 
   bool        colFirst, 
   bool        doBoth) 
   : m_name(name)
   , m_colFirst(colFirst)
   , m_doBoth(doBoth)
{
   assert(SPxScaler::isConsistent());
}

SPxScaler::SPxScaler(const SPxScaler& old)
   : m_name(old.m_name)
   , m_colscale(old.m_colscale)
   , m_rowscale(old.m_rowscale)
   , m_colFirst(old.m_colFirst)
   , m_doBoth(old.m_doBoth)
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

void SPxScaler::scale(SPxLP& lp) 
{
   int i;

   setup(lp);

   if (m_colFirst)
   {
      for(i = 0; i < lp.nCols(); ++i )
      {
         SVector& vec = lp.colVector_w(i);
         Real     x   = computeColscale(vec); 

         if (isZero(x))
            m_colscale[i] = 1.0;
         else
         {
            Real y          = 1.0 / x;
            m_colscale[i]   = y;
            vec            *= y;
            lp.maxObj_w(i) *= y;

            if (lp.upper(i) < infinity)
               lp.upper_w(i) *= x;
            if (lp.lower(i) > -infinity)
               lp.lower_w(i) *= x;
         }
      }
      
      for(i = 0; i < lp.nRows(); ++i )
      {
         SVector& vec = lp.rowVector_w(i);

         for(int j = 0; j < vec.size(); ++j )
            vec.value(j) *= m_colscale[vec.index(j)];

         Real x = computeRowscale(vec);

         if (isZero(x) || !m_doBoth)
            m_rowscale[i] = 1.0;
         else
         {
            Real y         = 1.0 / x;
            m_rowscale[i]  = y;
            vec           *= y;
            
            if (lp.rhs(i) < infinity)
               lp.rhs_w(i) *= y;
            if (lp.lhs(i) > -infinity)
               lp.lhs_w(i) *= y;
         }
      }
      if (m_doBoth)
      {
         for(i = 0; i < lp.nCols(); ++i )
         {
            SVector& vec = lp.colVector_w(i);
            
            for(int j = 0; j < vec.size(); ++j)
               vec.value(j) *= m_rowscale[vec.index(j)];
         }
      }
   }
   else
   {
      for(i = 0; i < lp.nRows(); ++i )
      {
         SVector& vec = lp.rowVector_w(i);
         Real     x   = computeRowscale(vec); // vec.maxAbs();

         if (isZero(x))
            m_rowscale[i] = 1.0;
         else
         {
            Real y        = 1.0 / x;
            m_rowscale[i] = y;
            vec          *= y;

            if (lp.rhs(i) < infinity)
               lp.rhs_w(i) *= y;
            if (lp.lhs(i) > -infinity)
               lp.lhs_w(i) *= y;
         }
      }
      for(i = 0; i < lp.nCols(); ++i )
      {
         SVector& vec = lp.colVector_w(i);

         for( int j = 0; j < vec.size(); ++j)
            vec.value(j) *= m_rowscale[vec.index(j)];

         Real x = computeColscale(vec);

         if (isZero(x) || !m_doBoth)
            m_colscale[i] = 1.0;
         else
         {
            Real y          = 1.0 / x;
            m_colscale[i]   = y;
            vec            *= y;
            lp.maxObj_w(i) *= y;

            if (lp.upper(i) < infinity)
               lp.upper_w(i) *= x;
            if (lp.lower(i) > -infinity)
               lp.lower_w(i) *= x;
         }
      }
      if (m_doBoth)
      {
         for( i = 0; i < lp.nRows(); ++i )
         {
            SVector& vec = lp.rowVector_w(i);
            
            for( int j = 0; j < vec.size(); ++j)
               vec.value(j) *= m_colscale[vec.index(j)];
         }
      }
   }
   assert(lp.isConsistent());
}

void SPxScaler::unscalePrimal(Vector& x) const
{
   assert(x.dim() == m_colscale.size());

   for(int i = 0; i < x.dim(); ++i )
      x[i] *= m_colscale[i];
}

void SPxScaler::unscaleDual(Vector& pi) const
{
   assert(pi.dim() == m_rowscale.size());

   ///@todo is this correct ?
   for(int i = 0; i < pi.dim(); ++i )
      pi[i] *= m_rowscale[i];
}

Real SPxScaler::minAbsColscale() const
{
   Real mini = infinity;

   for(int i = 0; i < m_colscale.size(); ++i)
      if (fabs(m_colscale[i]) < mini)
         mini = fabs(m_colscale[i]);

   return mini;
}

Real SPxScaler::maxAbsColscale() const
{
   Real maxi = 0.0;

   for(int i = 0; i < m_colscale.size(); ++i)
      if (fabs(m_colscale[i]) > maxi)
         maxi = fabs(m_colscale[i]);

   return maxi;
}

Real SPxScaler::minAbsRowscale() const
{
   Real mini = infinity;

   for(int i = 0; i < m_rowscale.size(); ++i)
      if (fabs(m_rowscale[i]) < mini)
         mini = fabs(m_rowscale[i]);

   return mini;
}

Real SPxScaler::maxAbsRowscale() const
{
   Real maxi = 0.0;

   for(int i = 0; i < m_rowscale.size(); ++i)
      if (fabs(m_rowscale[i]) > maxi)
         maxi = fabs(m_rowscale[i]);

   return maxi;
}

bool SPxScaler::isConsistent() const
{
   return m_colscale.isConsistent() && m_rowscale.isConsistent();
}

} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------


