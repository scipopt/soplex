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
#pragma ident "@(#) $Id: spxequilisc.cpp,v 1.9 2005/07/13 19:05:32 bzforlow Exp $"

/**@file  spxequilisc.cpp
 * @brief Equilibrium row/column scaling.
 */
#include <assert.h>

#include "spxequilisc.h"
#include "spxout.h"

namespace soplex
{
static const char* makename(bool colFirst, bool doBoth)
{
   const char* name;

   if (doBoth)
      name = colFirst ? "CR-Equilibrium" : "RC-Equilibrium";
   else
      name = colFirst ? "C-Equilibrium" : "R-Equilibrium";

   return name;
}

SPxEquiliSC::SPxEquiliSC(bool colFirst, bool doBoth)
   : SPxScaler(makename(colFirst, doBoth), colFirst, doBoth)
{}

Real SPxEquiliSC::computeScale(Real /*mini*/, Real maxi) const
{
   METHOD( "SPxEquiliSC::computeScale()" );

   return maxi;
}

void SPxEquiliSC::scale(SPxLP& lp) 
{
   METHOD( "SPxEquiliSC::scale()" );

   VERBOSE1({ s_spxout << "IEQUSC01 Equilibrium scaling LP" << std::endl; });   

   setup(lp);

   /* We want to do that direction first, with the lower ratio.
    * Reason:           
    *                               Rowratio
    *            0.04  0.02  0.01      4
    *            4000    20  1000    200
    * Colratio    1e5   1e3   1e5
    *
    * Row first =>                  Col next =>
    *               1   0.5  0.25         1   1   1 
    *               1   0.05 0.25         1  0.1  1
    *
    * Col first =>                  Row next =>
    *            1e-5  1e-3  1e-5        0.01  1  0.01
    *               1     1     1          1   1    1
    *
    */
   Real colratio = maxColRatio(lp);
   Real rowratio = maxRowRatio(lp);

   m_colFirst = colratio < rowratio;

   VERBOSE2({ s_spxout << "IEQUSC02 LP scaling statistics:" 
                       << " min= " << lp.minAbsNzo()
                       << " max= " << lp.maxAbsNzo()
                       << " col-ratio= " << colratio 
                       << " row-ratio= " << rowratio
                       << std::endl; });
   if (m_colFirst)
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

   VERBOSE3({ s_spxout << "IEQUSC03 \tRow scaling min= " << minAbsRowscale()
                       << " max= " << maxAbsRowscale()
                       << std::endl
                       << "\tCol scaling min= " << minAbsColscale()
                       << " max= " << maxAbsColscale()
                       << std::endl; });

   VERBOSE2({ s_spxout << "IEQUSC04 LP scaling statistics:" 
                       << " min= " << lp.minAbsNzo()
                       << " max= " << lp.maxAbsNzo()
                       << " col-ratio= " << maxColRatio(lp) 
                       << " row-ratio= " << maxRowRatio(lp) 
                       << std::endl; });
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





