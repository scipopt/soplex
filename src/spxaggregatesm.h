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
#pragma ident "@(#) $Id: spxaggregatesm.h,v 1.3 2001/11/22 08:57:22 bzfkocht Exp $"

#ifndef _SPXAGGREGATESM_H_
#define _SPXAGGREGATESM_H_

#include <assert.h>

#include "spxsimplifier.h"

namespace soplex
{
/** Remove redundant row and columns.
    !!! Das ist vermutlich falsch.
    This \Ref{SPxSimplifier} thries to eliminat redundant rows or columns from
    its loaded \Ref{SPxLP}.
 */
class SPxAggregateSM : public SPxSimplifier
{
protected:
   int eliminate(const SVector& row, double b);

public:
   ///
   double maxFill;
   ///
   double stability;
   ///
   void load(SPxLP*);
   ///
   int simplify();
   ///
   void unsimplify();
   ///
   double value(double x)
   {
      return x + lp->spxSense()*delta;
   }
   int isConsistent() const
   {
      return 1;
   };

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
