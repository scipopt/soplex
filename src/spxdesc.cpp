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
#pragma ident "@(#) $Id: spxdesc.cpp,v 1.1 2001/11/06 16:18:32 bzfkocht Exp $"

#include <iostream>
#include "spxbasis.h"

namespace soplex
{

//@ -----------------------------------------------------------------------------
/*      \Section{Complex Method Implementations}
        \SubSection{Descriptor}
 */

void SPxBasis::Desc::reSize(int rowDim, int colDim)
{
   rowstat.reSize(rowDim);
   colstat.reSize(colDim);
}

int SPxBasis::Desc::isConsistent() const
{
   return rowstat.isConsistent() && colstat.isConsistent();
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
