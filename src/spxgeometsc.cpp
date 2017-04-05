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

namespace soplex
{
/**@param maxIters   arbitrary small number, we choose 8
   @param minImpr    Bixby said Fourer said in MP 23, 274 ff. that 0.9 is a good value.
   @param goodEnough if the max/min ratio is already less then 1000/1 we do not scale.
*/ 

static Real computeScalingVec(
      const SVSet*           vecset,
      const DataArray<Real>& coScaleval,
      DataArray<Real>&       scaleval)
   {

      Real pmax = 0.0;

      for( int i = 0; i < vecset->num(); ++i )
      {
         const SVector& vec = (*vecset)[i];

         Real maxi = 0.0;
         Real mini = infinity;

         for( int j = 0; j < vec.size(); ++j )
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

         scaleval[i] = 1.0 / spxSqrt(mini * maxi);

         Real p = maxi / mini;

         if (p > pmax)
            pmax = p;
      }
      return pmax;
   }


SPxGeometSC::SPxGeometSC(int maxIters, Real minImpr, Real goodEnough)
   : SPxScaler("Geometric")
   , m_maxIterations(maxIters)
   , m_minImprovement(minImpr)
   , m_goodEnoughRatio(goodEnough)
{}

SPxGeometSC::SPxGeometSC(const SPxGeometSC& old)
   : SPxScaler(old)
   , m_maxIterations(old.m_maxIterations)
   , m_minImprovement(old.m_minImprovement)
   , m_goodEnoughRatio(old.m_goodEnoughRatio)
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

   Real pstart = 0.0;
   Real p0     = 0.0;
   Real p1     = 0.0;

   setup(lp);

   /* We want to do that direction first, with the lower ratio.
    * See SPxEquiliSC::scale() for a reasoning.
    */
   Real colratio = maxColRatio(lp);
   Real rowratio = maxRowRatio(lp);

   DataArray < Real >  rowscale(lp.nRows(), lp.nRows(), 1.2);
   DataArray < Real >  colscale(lp.nCols(), lp.nCols(), 1.2);

   bool colFirst = colratio < rowratio;

   MSG_INFO2( (*spxout), (*spxout) << "before scaling:"
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << colratio
                        << " row-ratio= " << rowratio
                        << std::endl; )

   // We make at most maxIterations.
   for( int count = 0; count < m_maxIterations; count++ )
   {
      if (colFirst)
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

      // record start value, this is done with m_col/rowscale = 1.0, so it is the
      // value from the "original" (as passed to the scaler) LP.
      if (count == 0)
      {
         pstart = p0;
         // are we already good enough ?
         if (pstart < m_goodEnoughRatio)
            break;
      }
      else // do not test at the first iteration, then abort if no improvement.
         if (p1 > m_minImprovement * p0)
            break;
   }

   // we scale only if either:
   // - we had at the beginning a ratio worse than 1000/1
   // - we have at least a 15% improvement.
   if( pstart < m_goodEnoughRatio || p1 > pstart * m_minImprovement )
   {
      MSG_INFO2( (*spxout), (*spxout) << "No scaling done." << std::endl; )
   }
   else
   {
      DataArray < int > colscaleExp = *m_activeColscaleExp;
      DataArray < int > rowscaleExp = *m_activeRowscaleExp;

      int i;
      for( i = 0; i < lp.nCols(); ++i )
      {
          frexp(colscale[i], &(colscaleExp[i]));
          colscaleExp[i] -= 1;
      }

      for( i = 0; i < lp.nRows(); ++i )
      {
          frexp(rowscale[i], &(rowscaleExp[i]));
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
