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
#pragma ident "@(#) $Id: spxequilisc.cpp,v 1.5 2003/01/10 12:46:14 bzfkocht Exp $"

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


Real SPxEquiliSC::computeColscale(const SVector& col) const
{
   return col.maxAbs();
}

Real SPxEquiliSC::computeRowscale(const SVector& row) const
{
   return row.maxAbs();
}

void SPxEquiliSC::scale(SPxLP& lp) 
{
   VERBOSE2({ std::cout << "Equilibrium scaling LP" << std::endl; });   

   VERBOSE3({ std::cout << "LP scaling statistics:" 
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << maxColRatio(lp) 
                        << " row-ratio= " << maxRowRatio(lp) 
                        << std::endl; });

   SPxScaler::scale(lp);

   VERBOSE3({ std::cout << "\tRow scaling min= " << minAbsRowscale()
                        << " max= " << maxAbsRowscale()
                        << std::endl
                        << "\tCol scaling min= " << minAbsColscale()
                        << " max= " << maxAbsColscale()
                        << std::endl; });

   VERBOSE3({ std::cout << "LP scaling statistics:" 
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << maxColRatio(lp) 
                        << " row-ratio= " << maxRowRatio(lp) 
                        << std::endl; });
}

#if 0
void SPxEquiliSC::scale(SPxLP& lp) 
{
   VERBOSE2({ std::cout << "Equilibrium scaling LP" << std::endl; });   

   VERBOSE3({ std::cout << "\tNon zero min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << std::endl; });
   int i;

   setup(lp);

   if (m_colFirst)
   {
      for(i = 0; i < lp.nCols(); ++i )
      {
         SVector& vec = lp.colVector_w(i);
         Real     x   = vec.maxAbs();

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
         Real     x   = 0.0;
         Real     y;

         for(int j = 0; j < vec.size(); ++j )
         {
            vec.value(j) *= m_colscale[vec.index(j)];
            y             = fabs(vec.value(j));
            x             = (x < y) ? y : x;
         }
#if 0
         if (lp.rhs(i) < infinity)
         {
            y = fabs(lp.rhs(i));
            x = (x < y) ? y : x;
         }
         if (lp.lhs(i) > -infinity)
         {
            y = fabs(lp.lhs(i));
            x = (x < y) ? y : x;
         }
#endif
         if (isZero(x) || !m_doBoth)
            m_rowscale[i] = 1.0;
         else
         {
            y              = 1.0 / x;
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
         Real     x   = vec.maxAbs();
         Real     y;
#if 0
         if (lp.rhs(i) < infinity)
         {
            y = fabs(lp.rhs(i));
            x = (x < y) ? y : x;
         }
         if (lp.lhs(i) > -infinity)
         {
            y = fabs(lp.lhs(i));
            x = (x < y) ? y : x;
         }
#endif
         if (isZero(x))
            m_rowscale[i] = 1.0;
         else
         {
            y             = 1.0 / x;
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
         Real     x   = 0;
         Real     y;

         for( int j = 0; j < vec.size(); ++j)
         {
            vec.value(j) *= m_rowscale[vec.index(j)];
            y             = fabs(vec.value(j));
            x             = (x < y) ? y : x;
         }
         if (isZero(x) || !m_doBoth)
            m_colscale[i] = 1.0;
         else
         {
            y               = 1.0 / x;
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



