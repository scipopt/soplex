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
#pragma ident "@(#) $Id: spxscaler.cpp,v 1.1 2002/04/04 14:59:04 bzfkocht Exp $"

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
   int i;

   s << sc.getName() << " scaler:";

   if (sc.m_lp == 0)
      s << "not initialised" << std::endl;
   else
   {
      s << std::endl;
      s << "colscale = [ ";
      for( i = 0; i < sc.m_lp->nCols(); i++ )
         s << sc.m_colscale[i] << " ";
      s << "]" << std::endl;
      
      s << "rowscale = [ ";
      for( i = 0; i < sc.m_lp->nRows(); i++ )
         s << sc.m_rowscale[i] << " ";
      s << "]" << std::endl;
   }
   return s;
}

SPxScaler::SPxScaler(
   const char* name, 
   bool        colFirst, 
   bool        doBoth) 
   : m_name(name)
   , m_lp(0)
   , m_colFirst(colFirst)
   , m_doBoth(doBoth)
{}

SPxScaler::~SPxScaler()
{
   m_name = 0;
   m_lp   = 0;
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

void SPxScaler::setLP(SPxLP* lp)
{
   assert(lp != 0);

   if (m_lp != lp)
   {
      m_lp = lp;

      m_colscale.reSize(m_lp->nCols());
      m_rowscale.reSize(m_lp->nRows());
   }
   int i;

   for( i = 0; i < m_lp->nCols(); i++ )
      m_colscale[i] = 1.0;

   for( i = 0; i < m_lp->nRows(); i++ )
      m_rowscale[i] = 1.0;
}

void SPxScaler::unscale()
{
   assert(m_lp != 0);
   assert(m_lp->isConsistent());

   int i;
   int j;

   for( i = 0; i < m_lp->nRows(); ++i )
   {
      SVector& vec = m_lp->rowVector_w(i);

      for( j = 0; j < vec.size(); ++j )
         vec.value(j) /= m_colscale[vec.index(j)];

      double y = 1.0 / m_rowscale[i];

      vec *= y;

      if (m_lp->rhs(i) < infinity)
         m_lp->rhs(i) *= y;

      if (m_lp->lhs(i) > -infinity)
         m_lp->lhs(i) *= y;
   }
   for( i = 0; i < m_lp->nCols(); ++i )
   {
      SVector& vec = m_lp->colVector_w(i);

      vec *= 1.0 / m_colscale[i];

      for( j = 0; j < vec.size(); ++j )
         vec.value(j) /= m_rowscale[vec.index(j)];

      m_lp->maxObj(i) /= m_colscale[i];

      if (m_lp->upper(i) < infinity)
         m_lp->upper(i) *= m_colscale[i];

      if (m_lp->lower(i) > -infinity)
         m_lp->lower(i) *= m_colscale[i];
   }
   assert(m_lp->isConsistent());
}

void SPxScaler::unscaleColVector(Vector& vec) const
{
   assert(m_lp      != 0);
   assert(vec.dim() == m_lp->nCols());

   for( int i = 0; i < m_lp->nCols(); ++i )
      vec[i] *= m_colscale[i];
}

void SPxScaler::unscaleColVector(SVector& vec) const
{
   assert(m_lp != 0);

   for( int i = 0; i < vec.size(); ++i )
   {
      assert(vec.index(i) < m_lp->nCols());

      vec.value(i) *= m_colscale[vec.index(i)];
   }
}

void SPxScaler::unscaleRowVector(Vector& vec) const
{
   assert(m_lp      != 0);
   assert(vec.dim() == m_lp->nRows());

   for( int i = 0; i < m_lp->nRows(); ++i )
      vec[i] *= m_rowscale[i];
}

void SPxScaler::unscaleRowVector(SVector& vec) const
{
   assert(m_lp != 0);

   for( int i = 0; i < vec.size(); ++i )
   {
      assert(vec.index(i) < m_lp->nRows());

      vec.value(i) *= m_rowscale[vec.index(i)];
   }
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


