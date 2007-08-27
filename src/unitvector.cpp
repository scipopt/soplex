/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: unitvector.cpp,v 1.8 2007/08/27 15:35:14 bzfberth Exp $"

#include <assert.h>

#include "unitvector.h"
#include "message.h"

namespace soplex
{

#ifndef NO_CONSISTENCY_CHECKS
bool UnitVector::isConsistent() const
{
   if (mem() != &themem)
      return MSGinconsistent("UnitVector");
   if (mem() + 1 != &themem1)
      return MSGinconsistent("UnitVector");
   if (size() != 1)
      return MSGinconsistent("UnitVector");
   if (max() != 1)
      return MSGinconsistent("UnitVector");

   return SVector::isConsistent();
}
#endif
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
