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
#pragma ident "@(#) $Id: spxrem1sm.h,v 1.5 2002/01/19 18:59:17 bzfkocht Exp $"

/**@file  spxrem1sm.h
 * @brief Remove singletons from LP.
 */
#ifndef _SPXREM1SM_H_
#define _SPXREM1SM_H_

#include <assert.h>

#include "real.h"
#include "spxsimplifier.h"

namespace soplex
{
/**@brief   LP simplifier for removing row/column singletons.
   @ingroup Algo

   This #SPxSimplifier removes rows and possibly columns containing one
   nonzero value only. Also empty rows and columns are removed.
*/
class SPxRem1SM : public SPxSimplifier
{
public:
   /// Remove singletons from the LP.
   int simplify();

   /// Reverse the doings of #simplify().
   void unsimplify();

   /// objective value for unsimplified LP.
   /**@todo This is different implementet then the others. 
    *       Why? Is this ok?
    */
   Real value(Real x)
   {
      return x - lp->spxSense() * delta;
   }
};
} // namespace soplex
#endif // _SPXREM1SM_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
