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
#pragma ident "@(#) $Id: unitvector.cpp,v 1.1 2001/11/06 16:18:32 bzfkocht Exp $"

/*      \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <assert.h>
#include <iostream>


/*  and class header files
 */
#include "unitvector.h"

namespace soplex
{


//@ ----------------------------------------------------------------------------
/*      \SubSection{Consistency}
 */
#define inconsistent                                                    \
{                                                                       \
std::cout << "ERROR: Inconsistency detected in class UnitVector\n";      \
return 0;                                                          \
}

int UnitVector::isConsistent() const
{
   if (mem() != &themem)
      inconsistent;
   if (mem() + 1 != &themem1)
      inconsistent;
   if (size() != 1)
      inconsistent;
   if (max() != 1)
      inconsistent;
   return SVector::isConsistent();
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
