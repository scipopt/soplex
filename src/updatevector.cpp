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
#pragma ident "@(#) $Id: updatevector.cpp,v 1.1 2001/11/06 16:18:32 bzfkocht Exp $"


/* \Section{Complex Methods}
 */

#include "updatevector.h"

namespace soplex
{

UpdateVector& UpdateVector::operator=(const UpdateVector& rhs)
{
   theval = rhs.theval;
   thedelta = rhs.thedelta;
   DVector::operator=(rhs);
   return *this;
}

#define inconsistent    \
{       \
std::cout << "Inconsistency detected in class UpdateVector\n";  \
return 0;  \
}

int UpdateVector::isConsistent() const
{
   if (dim() != thedelta.dim())
      inconsistent;
   return DVector::isConsistent()
          && thedelta.isConsistent();
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
