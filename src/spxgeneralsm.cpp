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
#pragma ident "@(#) $Id: spxgeneralsm.cpp,v 1.15 2002/04/14 12:41:54 bzfkocht Exp $"

//#define DEBUGGING 1

#include <iostream>

#include "spxdefines.h"
#include "spxgeneralsm.h"

namespace soplex
{
void SPxGeneralSM::load(SPxLP* p_lp)
{
   lp         = p_lp;
   m_rem1.load (p_lp);
   m_redu.load (p_lp);
   m_aggr.load (p_lp);
}

void SPxGeneralSM::unload()
{
   m_rem1.unload ();
   m_redu.unload ();
   m_aggr.unload ();
}

int SPxGeneralSM::simplify()
{
   int ret;
   int rows = lp->nRows();
   int cols = lp->nCols();
   int cnt  = rows + cols;

   if ((ret = m_rem1.simplify()) != 0) 
      return ret;
   if ((ret = m_redu.simplify()) != 0) 
      return ret;

   while(m_repth * cnt > lp->nRows() + lp->nCols())
   {
      cnt = lp->nRows() + lp->nCols();

      if ((ret = m_rem1.simplify()) != 0) 
         return ret;

      if (cnt == lp->nRows() + lp->nCols())
         break;

      // if ((ret = m_aggr.simplify()) != 0) 
      //   return ret;

      if ((ret = m_redu.simplify()) != 0) 
         return ret;
   }
   assert(lp->isConsistent());

   return 0;
}

/**@todo This is not correctly implented, since the simplifiers may be
 *      called several times one after the others, this sequence has
 *      to be tracked to make the calls of unsimplify in exactly the
 *      reverse order.
 *      At the moment this is irrelevant because all the unsimplifiers
 *      are not implemented anyway.
 */
void SPxGeneralSM::unsimplify()
{
   //m_rem1.unsimplify ();
   //m_aggr.unsimplify ();
   //m_redu.unsimplify ();
}

Real SPxGeneralSM::value(Real x)
{
   return m_rem1.value(m_aggr.value(m_redu.value(x)));
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
