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
#pragma ident "@(#) $Id: spxrem1sm.h,v 1.3 2001/11/22 08:57:24 bzfkocht Exp $"


#ifndef _SPXREM1SM_H_
#define _SPXREM1SM_H_

#include <assert.h>

#include "spxsimplifier.h"

namespace soplex
{
/** LP simplifier for removing singletons.
    This #SPxSimplifier removes rows and possibly columns containing one
    nonzero value only.
 */
class SPxRem1SM : public SPxSimplifier
{
public:
   ///
   int simplify();
   ///
   void unsimplify();
   ///
   double value(double x)
   {
      /**@todo Stimmt das und warum ist es anders als bei den anderen? */
      return x - lp->spxSense()*delta;
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
