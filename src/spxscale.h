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
#pragma ident "@(#) $Id: spxscale.h,v 1.3 2001/11/22 08:57:24 bzfkocht Exp $"

#ifndef _SPXSCALE_H_
#define _SPXSCALE_H_

#include <assert.h>

#include "spxsimplifier.h"
#include "dataarray.h"

namespace soplex
{
/** LP scaling.
    This #SPxSimplifier implementation performs simple scaling of the LPs rows
    and columns.
 */
class SPxScale : public SPxSimplifier
{
   DataArray < double > colscale;
   DataArray < double > rowscale;
   int                  rowScale;
   
public:
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

   /// default constructor.
   SPxScale() 
      : rowScale(1)
   {}

   /// destructor.
   virtual ~SPxScale()
   {}
};
} // namespace soplex
#endif // _SPXSCALE_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
