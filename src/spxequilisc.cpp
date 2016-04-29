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

/**@file  spxequilisc.cpp
 * @brief Equilibrium row/column scaling.
 */
#include <assert.h>

#include "spxequilisc.h"
#include "spxout.h"
#include "spxlpbase.h"

namespace soplex
{
static const char* makename(bool doBoth)
{
   return doBoth ? "bi-Equilibrium" : "uni-Equilibrium";
}

SPxEquiliSC::SPxEquiliSC(bool doBoth)
   : SPxScaler(makename(doBoth), false, doBoth)
{}

SPxEquiliSC::SPxEquiliSC(const SPxEquiliSC& old)
   : SPxScaler(old)
{}

SPxEquiliSC& SPxEquiliSC::operator=(const SPxEquiliSC& rhs)
{
   if(this != &rhs)
   {
      SPxScaler::operator=(rhs);
   }

   return *this;
}

Real SPxEquiliSC::computeScale(Real /*mini*/, Real maxi) const
{

   return maxi;
}

void SPxEquiliSC::scale(SPxLPBase<Real>& lp)
{

   MSG_INFO1( (*spxout), (*spxout) << "Equilibrium scaling LP" << std::endl; )

   setup(lp);

   /* We want to do the direction first, which has a lower maximal ratio,
    * since the lowest value in the scaled matrix is bounded from below by
    * the inverse of the maximum ratio of the direction that is done first
    * Example:
    *                     Rowratio
    *            0.1  1   10
    *            10   1   10
    *
    * Colratio   100  1
    *
    * Row first =>         Col next =>
    *            0.1  1          0.1  1
    *            1    0.1        1    0.1
    *
    * Col first =>         Row next =>
    *            0.01 1          0.01 1
    *            1    1          1    1
    *
    */
   Real colratio = maxColRatio(lp);
   Real rowratio = maxRowRatio(lp);

   bool colFirst = colratio < rowratio;

   MSG_INFO2( (*spxout), (*spxout) << "LP scaling statistics:"
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << colratio 
                        << " row-ratio= " << rowratio
                        << std::endl; )
   if (colFirst)
   {
      computeScalingVecs(lp.colSet(), m_rowscale, m_colscale);

      if (m_doBoth)
         computeScalingVecs(lp.rowSet(), m_colscale, m_rowscale);
   }
   else
   {
      computeScalingVecs(lp.rowSet(), m_colscale, m_rowscale);

      if (m_doBoth)
         computeScalingVecs(lp.colSet(), m_rowscale, m_colscale);
   }
   applyScaling(lp);

   MSG_INFO3( (*spxout), (*spxout) << "Row scaling min= " << minAbsRowscale()
                        << " max= " << maxAbsRowscale()
                        << std::endl
                        << "\tCol scaling min= " << minAbsColscale()
                        << " max= " << maxAbsColscale()
                        << std::endl; )

   MSG_INFO2( (*spxout), (*spxout) << "LP scaling statistics:"
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << maxColRatio(lp) 
                        << " row-ratio= " << maxRowRatio(lp) 
                        << std::endl; )
}

} // namespace soplex
