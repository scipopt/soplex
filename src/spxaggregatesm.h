/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2001 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxaggregatesm.h,v 1.5 2001/11/22 16:30:01 bzfkocht Exp $"

/**@file  spxaggregatesm.h
 * @brief LP variable aggregation.
 */
#ifndef _SPXAGGREGATESM_H_
#define _SPXAGGREGATESM_H_

#include <assert.h>

#include "spxsimplifier.h"

namespace soplex
{
/** @brief
    @ingroup Algo

    This #SPxSimplifier does variable aggregation.
 */
class SPxAggregateSM : public SPxSimplifier
{
private:
   static const double stability  = 0.01;  ///< stability factor.
   
   double maxFill;                         ///< ???  

   /// ???
   int eliminate(const SVector& row, double b);

public:
   /// Aggregate variable.
   int simplify();

   /// Undo #simplify().
   void unsimplify();

   /// objective value for unsimplified LP.
   double value(double x)
   {
      return x + lp->spxSense()*delta;
   }
};
} // namespace soplex
#endif // _SPXAGGREGATESM_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
