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
#pragma ident "@(#) $Id: spxscale.h,v 1.1 2001/11/06 16:18:33 bzfkocht Exp $"

#ifndef _SPXSCALE_H_
#define _SPXSCALE_H_


//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */

#include "spxsimplifier.h"
#include "dataarray.h"

namespace soplex
{






//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** LP scaling.
    This #SPxSimplifier# implementation performs simple scaling of the LPs rows
    and columns.
 */
class SPxScale : public SPxSimplifier
{
   SPxLP* lp;
   DataArray < double > colscale;
   DataArray < double > rowscale;
   int rowScale;
   
public:
   ///
   virtual void load(SPxLP*);
   ///
   virtual void unload();
   ///
   virtual SPxLP* loadedLP() const
   {
      return lp;
   }
   ///
   virtual int simplify();
   ///
   virtual void unsimplify();

   ///
   int isConsistent() const
   {
      return colscale.isConsistent()
             && rowscale.isConsistent();
   }
   ///
   SPxScale() : rowScale(1)
   {}

   virtual ~SPxScale()
   {}

};


} // namespace soplex
#endif // _SPXSCALE_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
