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
#pragma ident "@(#) $Id: unitvector.cpp,v 1.3 2001/12/25 16:03:25 bzfkocht Exp $"

#include <assert.h>

#include "unitvector.h"
#include "spxmessage.h"

namespace soplex
{

int UnitVector::isConsistent() const
{
   if (mem() != &themem)
      return SPXinconsistent("UnitVector");
   if (mem() + 1 != &themem1)
      return SPXinconsistent("UnitVector");
   if (size() != 1)
      return SPXinconsistent("UnitVector");
   if (max() != 1)
      return SPXinconsistent("UnitVector");

   return SVector::isConsistent();
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
