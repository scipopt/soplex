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
#pragma ident "@(#) $Id: spxequilisc.cpp,v 1.6 2003/01/12 13:09:40 bzfkocht Exp $"

/**@file  spxequilisc.cpp
 * @brief Equilibrium row/column scaling.
 */
#include <assert.h>

#include "spxequilisc.h"

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

   VERBOSE1({ std::cout << "Equilibrium scaling LP" << std::endl; });   

   setup(lp);

   VERBOSE2({ std::cout << "LP scaling statistics:" 
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << maxColRatio(lp) 
                        << " row-ratio= " << maxRowRatio(lp) 
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

   VERBOSE3({ std::cout << "\tRow scaling min= " << minAbsRowscale()
                        << " max= " << maxAbsRowscale()
                        << std::endl
                        << "\tCol scaling min= " << minAbsColscale()
                        << " max= " << maxAbsColscale()
                        << std::endl; });

   VERBOSE2({ std::cout << "LP scaling statistics:" 
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



