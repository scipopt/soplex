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
#pragma ident "@(#) $Id: spxgeometsc.cpp,v 1.2 2003/01/10 12:46:14 bzfkocht Exp $"

/**@file  spxgeometsc.cpp
 * @brief Geometric mean row/column scaling.
 */
#include <assert.h>

#include "spxgeometsc.h"

namespace soplex
{
static const char* makename(bool colFirst, bool doBoth)
{
   const char* name;

   if (doBoth)
      name = colFirst ? "CR-Geometric" : "RC-Geometric";
   else
      name = colFirst ? "C-Geometric" : "R-Geometric";

   return name;
}

SPxGeometSC::SPxGeometSC(bool colFirst, bool doBoth)
   : SPxScaler(makename(colFirst, doBoth), colFirst, doBoth)
{}


Real SPxGeometSC::computeColscale(const SVector& col) const
{
   Real mini = col.minAbs();
   Real maxi = col.maxAbs();

   return sqrt(mini * maxi);
}

Real SPxGeometSC::computeRowscale(const SVector& row) const
{
   Real mini = row.minAbs();
   Real maxi = row.maxAbs();

   return sqrt(mini * maxi);
}


Real SPxGeometSC::doRow(const SPxLP& lp) 
{
   Real pmax = 0.0;

   for(int i = 0; i < lp.nRows(); ++i )
   {
      const SVector& vec = lp.rowVector(i);
            
      Real maxi = 0.0;
      Real mini = infinity;
            
      for( int j = 0; j < vec.size(); ++j)
      {
         Real x = fabs(vec.value(j) * m_colscale[vec.index(j)]);
               
         if (!isZero(x))
         {
            if (x > maxi)
               maxi = x;
            if (x < mini)
               mini = x;
         }
      }
      m_rowscale[i] = 1.0 / sqrt(mini * maxi);
            
      Real p = maxi / mini;
            
      if (p > pmax)
         pmax = p;
   }
   return pmax;
}

Real SPxGeometSC::doCol(const SPxLP& lp) 
{
   Real pmax = 0.0;

   for(int i = 0; i < lp.nCols(); ++i )
   {
      const SVector& vec = lp.colVector(i);
            
      Real maxi = 0.0;
      Real mini = infinity;
            
      for( int j = 0; j < vec.size(); ++j)
      {
         Real x = fabs(vec.value(j) * m_rowscale[vec.index(j)]);
         
         if (!isZero(x))
         {
            if (x > maxi)
               maxi = x;
            if (x < mini)
               mini = x;
         }
      }
      m_colscale[i] = 1.0 / sqrt(mini * maxi);
      
      Real p = maxi / mini;
      
      if (p > pmax)
         pmax = p;
   }
   return pmax;
}

#if 0
void SPxGeometSC::scale(SPxLP& lp) 
{
   VERBOSE2({ std::cout << "Geometric scaling LP" << std::endl; });   

   VERBOSE3({ std::cout << "\tNon zero min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << std::endl; });
   
   SPxScaler::scale(lp);

   VERBOSE3({ std::cout << "\tRow scaling min= " << minAbsRowscale()
                        << " max= " << maxAbsRowscale()
                        << std::endl
                        << "\tCol scaling min= " << minAbsColscale()
                        << " max= " << maxAbsColscale()
                        << std::endl
                        << "\tNon zero min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << std::endl; });
}
#else
void SPxGeometSC::scale(SPxLP& lp) 
{
   Real pstart = 0.0;
   Real p0;
   Real p1;

   setup(lp);

   VERBOSE2({ std::cout << "Geometric scaling LP" << std::endl; });   

   VERBOSE3({ std::cout << "LP scaling statistics:" 
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << maxColRatio(lp)
                        << " row-ratio= " << maxRowRatio(lp)
                        << std::endl; });

   // We make at most 8 (random number) iterations. 
   for(int count = 0; count < 8; count++)
   {
      if (m_colFirst)
      {
         p0 = doCol(lp);
         p1 = doRow(lp);
      }
      else
      {
         p0 = doRow(lp);
         p1 = doCol(lp);
      }
      VERBOSE3({ std::cout << "Geometric scaling round " << count
                           << " col-ratio= " << (m_colFirst ? p0 : p1)
                           << " row-ratio= " << (m_colFirst ? p1 : p0)
                           << std::endl; });

      // rcord start value
      if (count == 0)
         pstart = p0;
      else // do not test at the first iteration, then abort if no improvement.
         if (p1 > 0.9 * p0)
            break;
   }      
   
   // we scale only if either
   // we had at the beginng a ratio worse then 1/100
   // we have at least a 10% improvement.
   if (pstart < 1e3 || p1 > pstart * 0.9)
   {
      setup(lp);

      VERBOSE2({ std::cout << "No scaling done." << std::endl; });
   }
   else
   {
      int i;

      // now doit
      for(i = 0; i < lp.nRows(); ++i )
      {
         SVector& vec = lp.rowVector_w(i);

         for( int j = 0; j < vec.size(); ++j)
            vec.value(j) *= m_colscale[vec.index(j)] * m_rowscale[i];

         if (lp.rhs(i) < infinity)
            lp.rhs_w(i) *= m_rowscale[i];
         if (lp.lhs(i) > -infinity)
            lp.lhs_w(i) *= m_rowscale[i];
      }
      for(i = 0; i < lp.nCols(); ++i )
      {
         SVector& vec = lp.colVector_w(i);
         
         for( int j = 0; j < vec.size(); ++j)
            vec.value(j) *= m_rowscale[vec.index(j)] * m_colscale[i];
         
         lp.maxObj_w(i) *= m_colscale[i];
         
         if (lp.upper(i) < infinity)
            lp.upper_w(i) /= m_colscale[i];
         if (lp.lower(i) > -infinity)
            lp.lower_w(i) /= m_colscale[i];
      }
   }
   assert(lp.isConsistent());

   VERBOSE3({ std::cout << "Row scaling min= " << minAbsRowscale()
                        << " max= " << maxAbsRowscale()
                        << std::endl
                        << "Col scaling min= " << minAbsColscale()
                        << " max= " << maxAbsColscale()
                        << std::endl; });

   VERBOSE3({ std::cout << "LP scaling statistics:" 
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << maxColRatio(lp) 
                        << " row-ratio= " << maxRowRatio(lp) 
                        << std::endl; });
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



