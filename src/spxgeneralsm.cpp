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
#pragma ident "@(#) $Id: spxgeneralsm.cpp,v 1.17 2003/01/10 12:46:14 bzfkocht Exp $"

//#define DEBUGGING 1

#include <iostream>

#include "spxdefines.h"
#include "spxgeneralsm.h"

namespace soplex
{
SPxSimplifier::Result SPxGeneralSM::simplify(SPxLP& lp, Real eps, Real delta)
{
   Result ret;
#if 0
   int    rows = lp.nRows();
   int    cols = lp.nCols();
   int    cnt  = rows + cols;
#endif

   eps   *= 100.0;  // 1e-14 -> 1e-12
   delta *= 0.001;  // 1e-6  -> 1e-9

   if ((ret = m_inter.simplify(lp, eps, delta)) != OKAY)
      return ret;
   if ((ret = m_rem1.simplify(lp, eps, delta)) != OKAY) 
      return ret;
#if 0

   while(m_repth * cnt > lp.nRows() + lp.nCols())
   {
      cnt = lp.nRows() + lp.nCols();

      if ((ret = m_rem1.simplify()) != 0) 
         return ret;

      if (cnt == lp.nRows() + lp.nCols())
         break;

      // if ((ret = m_aggr.simplify()) != 0) 
      //   return ret;

      if ((ret = m_redu.simplify()) != 0) 
         return ret;
   }
#endif
   assert(lp.isConsistent());

   return OKAY;
}

const Vector& SPxGeneralSM::unsimplifiedPrimal(const Vector& x)
{
   return m_inter.unsimplifiedPrimal(m_rem1.unsimplifiedPrimal(x));
}

const Vector& SPxGeneralSM::unsimplifiedDual(const Vector& pi)
{
   return m_inter.unsimplifiedDual(m_rem1.unsimplifiedDual(pi));
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
