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
#pragma ident "@(#) $Id: spxrem1sm.h,v 1.1 2001/11/06 16:18:33 bzfkocht Exp $"


#ifndef _SPXREM1SM_H_
#define _SPXREM1SM_H_

//@ ----------------------------------------------------------------------------
/*  \Section{Imports}
    Import required system include files ...
 */
#include <assert.h>


/*  ... and class header files
 */

#include "spxsimplifier.h"

namespace soplex
{






//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** LP simplifier for removing singletons.
    This \Ref{SPxSimplifier} removes rows and possibly columns containing one
    nonzero value only.
 */
class SPxRem1SM : public SPxSimplifier
{
private:
   SPxLP* lp;
   double delta;

protected:
public:
   ///
   void load(SPxLP*);
   ///
   void unload();
   ///
   SPxLP* loadedLP() const
   {
      return lp;
   }
   ///
   int simplify();
   ///
   void unsimplify();
   ///
   double value(double x)
   {
      return x - lp->spxSense()*delta;
   }
};

} // namespace soplex
#endif // _SPXREM1SM_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
