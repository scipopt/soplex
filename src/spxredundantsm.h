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
#pragma ident "@(#) $Id: spxredundantsm.h,v 1.4 2001/11/22 16:30:01 bzfkocht Exp $"

/**@file  spxredundant.h
 * @brief Remove redundant row and columns.
 */
#ifndef _SPXREDUNDANTSM_H_
#define _SPXREDUNDANTSM_H_

#include <assert.h>

#include "spxsimplifier.h"

namespace soplex
{
/**@brief   Remove redundant row and columns.
   @ingroup Algo

   This #SPxSimplifier tries to eliminate redundant rows and columns from
   its loaded #SPxLP.
 */
class SPxRedundantSM : public SPxSimplifier
{
private:
   static const double eps = 1e-10; ///< epslion for what is regarded equal.

public:
   /// Remove redundant rows and columns.
   int simplify();

   /// Reverse the doings of #simplify().
   void unsimplify();

   /// objective value for unsimplified LP.
   double value(double x)
   {
      return x + lp->spxSense() * delta;
   }
};
} // namespace soplex
#endif // _SPXREDUNDANTSM_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
