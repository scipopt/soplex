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
#pragma ident "@(#) $Id: spxsimplifier.h,v 1.1 2001/11/06 16:18:33 bzfkocht Exp $"

#ifndef _SPXSIMPLIFIER_H_
#define _SPXSIMPLIFIER_H_


//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */

#include "spxlp.h"

namespace soplex
{






//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** LP simplification base class.
    Instances of classes derived from #SPxSimplifier# may be loaded to #SoPlex# in
    order to simplify LPs before solving them. #SoPlex# will #load()# itself to
    the #SPxSimplifier# and then call #simplify()#. Generally any #SPxLP# can be
    loaded to a #SPxSimplifier# for #simplify()#ing it. The simplification can
    be undone by calling #unsimplify()#.
 */
class SPxSimplifier
{
private:
protected:
public:
   ///
   virtual void load(SPxLP*) = 0;
   ///
   virtual void unload() = 0;
   ///
   virtual SPxLP* loadedLP() const = 0;

   /** Simplify loaded #SPxLP#. It returns
       \begin{description}
       \item[0]     if this could be done,
       \item[1]     if the LP was detected to be unbounded or
       \item[-1]    if the LP was detected to be infeasible.
       \end{description}
    */
   virtual int simplify() = 0;
   ///
   virtual void unsimplify() = 0;

   /** objective value for unsimplified LP.
       The simplifyed LP may show other objective values than the
       original, if a constant part has been removed from the LP.
       This method returns the value for the original LP, for a
       value #x# of the simplified LP.
    */
   virtual double value(double x)
   {
      return x;
   }

   ///
   int isConsistent() const
   {
      return 1;
   }
};


} // namespace soplex
#endif // _SPXSIMPLIFIER_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
