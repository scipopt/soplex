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
#pragma ident "@(#) $Id: unitvector.h,v 1.2 2001/11/06 23:31:07 bzfkocht Exp $"


#ifndef _UNITVECTOR_H_
#define _UNITVECTOR_H_

//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */

#include "svector.h"

namespace soplex
{






//@ ----------------------------------------------------------------------------
/*      \Section{Class Declaration}
 */

/** sparse unit vector.
    A #UnitVector# is an #SVector# that can take only one nonzero value with
    value 1 but arbitrary index.
 */
class UnitVector : public SVector
{
private:
   Element themem;
   Element themem1;

public:
   ///
   double value(int n) const
   {
      (void)n;
      assert(n == 0);
      return 1;
   }

   /// construct #i#-th unit vector.
   UnitVector(int i = 0)
      : SVector(2, &themem)
   {
      add(i, 1.0);
   }

   ///
   UnitVector(const UnitVector& rhs)
      : SVector(2, &themem)
         , themem (rhs.themem)
         , themem1(rhs.themem1)
   {}

   ///

   UnitVector& operator=(const UnitVector& rhs)
   {
      themem1 = rhs.themem1;
      return *this;
   }

   ///
   int isConsistent() const;
};


} // namespace soplex
#endif // _UNITVECTOR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
