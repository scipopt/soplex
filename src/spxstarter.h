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
#pragma ident "@(#) $Id: spxstarter.h,v 1.1 2001/11/06 16:18:33 bzfkocht Exp $"


#ifndef _SPXDSTARTER_H_
#define _SPXDSTARTER_H_

//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */

#include "soplex.h"

namespace soplex
{






//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** #SoPlex# start basis generation base class.
    #SPxStarter# is the virtual base class for classes generating a starter basis
    for the Simplex solver #SoPlex#. When a #SPxStarter# object has been loaded
    to a #SoPlex# solver, the latter will call method #generate()# in order to
    have a start basis generated. Implementations of method #generate()# must
    terminate by #SoPlex::load()#ing the generated basis to #SoPlex#. Loaded
    basises must be nonsingular.
 */
class SPxStarter
{
public:
   /// generate start basis for loaded basis.
   virtual void generate(SoPlex& base) = 0;

   /// destructor.
   virtual ~SPxStarter()
   { }

   /// check consistency.

   virtual int isConsistent() const;
};


} // namespace soplex
#endif // _SPXDSTARTER_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
