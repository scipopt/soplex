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
#pragma ident "@(#) $Id: spxscale.h,v 1.4 2001/11/22 16:30:01 bzfkocht Exp $"

/**@file  spxscale.h
 * @brief LP row/column scaling.
 */
#ifndef _SPXSCALE_H_
#define _SPXSCALE_H_

#include <assert.h>

#include "spxsimplifier.h"
#include "dataarray.h"

namespace soplex
{
/**@brief LP row/column scaling.
   @ingroup Algo

   This #SPxSimplifier implementation performs simple scaling of the 
   LPs rows and columns.

   @todo The type of scaling (row/column) can is hard coded. This should
         br selectable.
*/
class SPxScale : public SPxSimplifier
{
   DataArray < double > colscale;  ///< column scaleing factors
   DataArray < double > rowscale;  ///< row scaleing factors
   bool                 rowScale;  ///< do row scaleing (not column scaleing)
   
public:
   /// Scale the loaded #SPxLP.
   virtual int simplify();

   /// Unscale the loaded #SPxLP.
   virtual void unsimplify();

   /// consistency check.
   int isConsistent() const
   {
      return colscale.isConsistent() && rowscale.isConsistent();
   }

   /// default constructor.
   SPxScale() 
      : rowScale(true)
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
