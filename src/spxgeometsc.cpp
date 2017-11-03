/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  spxgeometsc.cpp
 * @brief Geometric mean row/column scaling.
 */
#include <assert.h>

#include "spxgeometsc.h"
#include "spxout.h"
#include "spxlpbase.h"
#include <vector>

namespace soplex
{

static Real computeScalingVec(
      const SVSet*             vecset,
      const std::vector<Real>& coScaleval,
      std::vector<Real>&       scaleval)
   {

      Real pmax = 0.0;
      for( int i = 0; i < vecset->num(); ++i )
      {
         const SVector& vec = (*vecset)[i];

         Real maxi = 0.0;
         Real mini = infinity;

         for( int j = 0; j < vec.size(); ++j )
         {
            const Real x = spxAbs(vec.value(j) * coScaleval[unsigned(vec.index(j))]);

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

         scaleval[unsigned(i)] = 1.0 / spxSqrt(mini * maxi);

         const Real p = maxi / mini;

         if (p > pmax)
            pmax = p;
      }
      return pmax;
   }


SPxGeometSC::SPxGeometSC(int maxIters, Real minImpr, Real goodEnough, bool equilibrate)
   : SPxScaler("Geometric")
   , m_maxIterations(maxIters)
   , m_minImprovement(minImpr)
   , m_goodEnoughRatio(goodEnough)
   , postequilibration(equilibrate)
{}

SPxGeometSC::SPxGeometSC(const SPxGeometSC& old)
   : SPxScaler(old)
   , m_maxIterations(old.m_maxIterations)
   , m_minImprovement(old.m_minImprovement)
   , m_goodEnoughRatio(old.m_goodEnoughRatio)
   , postequilibration(old.postequilibration)
{}

SPxGeometSC& SPxGeometSC::operator=(const SPxGeometSC& rhs)
{
   if (this != &rhs)
   {
      SPxScaler::operator=(rhs);
   }

   return *this;
}

void SPxGeometSC::scale(SPxLPBase<Real>& lp, bool persistent)
{

   MSG_INFO1( (*spxout), (*spxout) << "Geometric scaling LP" << (persistent ? " (persistent)" : "") << std::endl; )

   setup(lp);

   /* We want to do that direction first, with the lower ratio.
    * See SPxEquiliSC::scale() for a reasoning.
    */
   const Real colratio = maxColRatio(lp);
   const Real rowratio = maxRowRatio(lp);

   bool colFirst = colratio < rowratio;

   Real p0start;
   Real p1start;

   if( colFirst )
   {
     p0start = colratio;
     p1start = rowratio;
   }
   else
   {
     p0start = rowratio;
     p1start = colratio;
   }

   MSG_INFO2( (*spxout), (*spxout) << "before scaling:"
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << colratio
                        << " row-ratio= " << rowratio
                        << std::endl; )

   // are we already good enough ?
   if( p1start < m_goodEnoughRatio )
   {
      MSG_INFO2( (*spxout), (*spxout) << "No scaling done." << std::endl; )
      lp.setScalingInfo(true);
      return;
   }

   std::vector<Real> rowscale(unsigned(lp.nRows()), 1.0);
   std::vector<Real> colscale(unsigned(lp.nCols()), 1.0);

   Real p0 = 0.0;
   Real p1 = 0.0;
   Real p0prev = p0start;
   Real p1prev = p1start;

   // We make at most maxIterations.
   for( int count = 0; count < m_maxIterations; count++ )
   {
      if( colFirst )
      {
         p0 = computeScalingVec(lp.colSet(), rowscale, colscale);
         p1 = computeScalingVec(lp.rowSet(), colscale, rowscale);
      }
      else
      {
         p0 = computeScalingVec(lp.rowSet(), colscale, rowscale);
         p1 = computeScalingVec(lp.colSet(), rowscale, colscale);
      }

      MSG_INFO3( (*spxout), (*spxout) << "Geometric scaling round " << count
                           << " col-ratio= " << (colFirst ? p0 : p1)
                           << " row-ratio= " << (colFirst ? p1 : p0)
                           << std::endl; )

      if( p0 > m_minImprovement * p0prev && p1 > m_minImprovement * p1prev )
            break;

      p0prev = p0;
      p1prev = p1;
   }

   // we scale only if we have enough (15%) improvement.
   if( p0 > m_minImprovement * p0start && p1 > m_minImprovement * p1start )
   {
      MSG_INFO2( (*spxout), (*spxout) << "No scaling done." << std::endl; )
      lp.setScalingInfo(true);
   }
   else
   {
      DataArray < int >& colscaleExp = *m_activeColscaleExp;
      DataArray < int >& rowscaleExp = *m_activeRowscaleExp;

      for( int i = 0; i < lp.nCols(); ++i )
      {
          frexp(double(colscale[unsigned(i)]), &(colscaleExp[i]));
          colscaleExp[i] -= 1;
      }

      for( int i = 0; i < lp.nRows(); ++i )
      {
          frexp(double(rowscale[unsigned(i)]), &(rowscaleExp[i]));
          rowscaleExp[i] -= 1;
      }

      applyScaling(lp);

      MSG_INFO3( (*spxout), (*spxout) << "Row scaling min= " << minAbsRowscale()
                           << " max= " << maxAbsRowscale()
                           << std::endl
                           << "IGEOSC06 Col scaling min= " << minAbsColscale()
                           << " max= " << maxAbsColscale()
                           << std::endl; )

      MSG_INFO2( (*spxout), (*spxout) << "after scaling: "
                           << " min= " << lp.minAbsNzo(false)
                           << " max= " << lp.maxAbsNzo(false)
                           << " col-ratio= " << maxColRatio(lp) 
                           << " row-ratio= " << maxRowRatio(lp) 
                           << std::endl; )
   }
}

} // namespace soplex
